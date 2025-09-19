from libra_toolbox.tritium.model import ureg, Model, quantity_to_activity
import numpy as np
import json
from libra_toolbox.tritium.lsc_measurements import (
    LIBRARun,
    LSCFileReader,
    GasStream,
    LSCSample,
    LIBRASample,
)
from pathlib import Path
from datetime import datetime


all_file_readers = []
all_quench = []


def create_sample(label: str, filename: str) -> LSCSample:
    """
    Create a LSCSample from a LSC file with background substracted.

    Args:
        label: the label of the sample in the LSC file
        filename: the filename of the LSC file

    Returns:
        the LSCSample object
    """
    # check if a LSCFileReader has been created for this filename
    found = False
    for file_reader in all_file_readers:
        if file_reader.filename == filename:
            found = True
            break

    # if not, create it and add it to the list of LSCFileReaders
    if not found:
        file_reader = LSCFileReader(filename, labels_column="SMPL_ID")

    file_reader.read_file()

    # create the sample
    sample = LSCSample.from_file(file_reader, label)

    # try to find the background sample from the file
    background_labels = ["1L-BL-1", "1L-BL-2", "1L-BL-3", "1L-BL-4"]
    background_sample = None

    for background_label in background_labels:
        try:
            background_sample = LSCSample.from_file(file_reader, background_label)
            break
        except ValueError:
            continue

    if background_sample is None:
        raise ValueError(f"Background sample not found in {filename}")

    # substract background
    sample.substract_background(background_sample)

    # read quench set
    all_quench.append(file_reader.quench_set)

    return sample


lsc_data_folder = "../../data/tritium_detection"
with open("../../data/general.json", "r") as f:
    general_data = json.load(f)

run_nb = general_data["general_data"]["run_nb"]


# read start time from general.json
all_start_times = []
for generator in general_data["generators"]:
    if generator["enabled"] is False:
        continue
    for irradiation_period in generator["periods"]:
        start_time = datetime.strptime(irradiation_period["start"], "%m/%d/%Y %H:%M")
        all_start_times.append(start_time)
start_time = min(all_start_times)


# create gas streams
gas_streams = {}
for stream, samples in general_data["tritium_detection"].items():
    stream_samples = []
    for sample_nb, sample_dict in samples.items():
        libra_samples = []
        if sample_dict["actual_sample_time"] is None:
            continue
        for vial_nb, filename in sample_dict["lsc_vials_filenames"].items():
            sample = create_sample(
                label=f"1L-{stream}_{run_nb}-{sample_nb}-{vial_nb}",
                filename=f"{lsc_data_folder}/{filename}",
            )
            libra_samples.append(sample)

        time_sample = datetime.strptime(
            sample_dict["actual_sample_time"], "%m/%d/%Y %H:%M"
        )
        stream_samples.append(LIBRASample(libra_samples, time=time_sample))
    gas_streams[stream] = GasStream(stream_samples, start_time=start_time)


# create run
run = LIBRARun(streams=list(gas_streams.values()), start_time=start_time)

# check that only one quench set is used
assert len(np.unique(all_quench)) == 1

# check that background is always substracted  # TODO this should be done automatically in LIBRARun
for stream in run.streams:
    for sample in stream.samples:
        for lsc_vial in sample.samples:
            assert (
                lsc_vial.background_substracted
            ), f"Background not substracted for {sample}"

IV_stream = gas_streams["IV"]
OV_stream = gas_streams["OV"]

sampling_times = {
    "IV": sorted(IV_stream.relative_times_as_pint),
    "OV": sorted(OV_stream.relative_times_as_pint),
}

replacement_times_top = sampling_times["IV"]
replacement_times_walls = sampling_times["OV"]

# read gas change time
if general_data["cover_gas"]["switched_to"]["gas_switch_time"]:
    gas_switch_time = datetime.strptime(
        general_data["cover_gas"]["switched_to"]["gas_switch_time"], "%m/%d/%Y %H:%M"
    )
    gas_switch_deltatime = gas_switch_time - start_time
    gas_switch_deltatime = gas_switch_deltatime.total_seconds() * ureg.s
    gas_switch_deltatime = gas_switch_deltatime.to(ureg.day)

# tritium model

baby_diameter = 14 * ureg.cm  # TODO confirm with CAD
baby_radius = 0.5 * baby_diameter
baby_volume = 1 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section

# read irradiation times from general.json

irradiations = []
for generator in general_data["generators"]:
    if generator["enabled"] is False:
        continue
    for irradiation_period in generator["periods"]:
        irr_start_time = (
            datetime.strptime(irradiation_period["start"], "%m/%d/%Y %H:%M")
            - start_time
        )
        irr_stop_time = (
            datetime.strptime(irradiation_period["end"], "%m/%d/%Y %H:%M") - start_time
        )
        irr_start_time = irr_start_time.total_seconds() * ureg.second
        irr_stop_time = irr_stop_time.total_seconds() * ureg.second
        irradiations.append([irr_start_time, irr_stop_time])

# Neutron rate


# check if neutron rate is provided in processed_data.json
processed_data_file = Path("../../data/processed_data.json")
neutron_rate = None
if processed_data_file.exists():
    with open(processed_data_file, "r") as f:
        processed_data = json.load(f)
    if "neutron_rate_used_in_model" in processed_data:
        if processed_data["neutron_rate_used_in_model"]["value"] is not None:
            neutron_rate = processed_data["neutron_rate_used_in_model"]["value"] * ureg(
                processed_data["neutron_rate_used_in_model"]["unit"]
            )
            neutron_rate_uncertainty = processed_data["neutron_rate_used_in_model"]["error"]* ureg(
                    processed_data["neutron_rate_used_in_model"]["unit"]
                )
            print(
                f"Using neutron rate from processed_data.json: {neutron_rate} ± {neutron_rate_uncertainty}"
            )
if neutron_rate is None:
    neutron_rate = 9.0e07 * ureg.neutron * ureg.s**-1 # based on manufacturer test data for generator settings
    neutron_rate_uncertainty = 2.0e07 * ureg.neutron * ureg.s**-1
    print(
        f"Using default neutron rate: {neutron_rate} ± {neutron_rate_uncertainty}"
    )

neutron_rate_relative_uncertainty = (neutron_rate_uncertainty / neutron_rate).to(
    ureg.dimensionless
)

# TBR from OpenMC

filename = "../neutron/statepoint.100.h5"
filename = Path(filename)

if not filename.exists():
    raise FileNotFoundError(f"{filename} does not exist, run OpenMC first")

import openmc

sp = openmc.StatePoint(filename)
tally_df = sp.get_tally(name="TBR").get_pandas_dataframe()
calculated_TBR = tally_df["mean"].iloc[0] * ureg.particle * ureg.neutron**-1
calculated_TBR_std_dev = (
    tally_df["std. dev."].iloc[0] * ureg.particle * ureg.neutron**-1
)

# TBR from measurements

total_irradiation_time = sum([irr[1] - irr[0] for irr in irradiations])

T_consumed = neutron_rate * total_irradiation_time

# to calculate the measured TBR we ignore the last samples for which
# we have some contribution from other sources (nGen, cyclotron, etc.)
nb_samples_included_iv = 9
T_produced_IV = IV_stream.get_cumulative_activity("total")[nb_samples_included_iv - 1]

measured_TBR = (T_produced_IV / quantity_to_activity(T_consumed)).to(
    ureg.particle * ureg.neutron**-1
)

# Run 1 transport coeff and measured TBR for overlay
optimised_ratio = 0
k_top = 12 * 1.45 * 8.9e-8 * ureg.m * ureg.s**-1
k_wall = optimised_ratio * k_top

baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=measured_TBR,
    neutron_rate=neutron_rate,
    irradiations=irradiations,
    k_top=k_top,
    k_wall=k_wall,
)


# store processed data
processed_data = {
    "modelled_baby_radius": {
        "value": baby_radius.magnitude,
        "unit": str(baby_radius.units),
    },
    "modelled_baby_height": {
        "value": baby_height.magnitude,
        "unit": str(baby_height.units),
    },
    "irradiations": [
        {
            "start_time": {
                "value": irr[0].magnitude,
                "unit": str(irr[0].units),
            },
            "stop_time": {
                "value": irr[1].magnitude,
                "unit": str(irr[1].units),
            },
        }
        for irr in irradiations
    ],
    "neutron_rate_used_in_model": {
        "value": baby_model.neutron_rate.magnitude,
        "unit": str(baby_model.neutron_rate.units),
        "error": neutron_rate_uncertainty.magnitude,
    },
    "measured_TBR": {
        "value": measured_TBR.magnitude,
        "unit": str(measured_TBR.units),
    },
    "TBR_used_in_model": {
        "value": baby_model.TBR.magnitude,
        "unit": str(baby_model.TBR.units),
    },
    "k_top": {
        "value": baby_model.k_top.magnitude,
        "unit": str(baby_model.k_top.units),
    },
    "k_wall": {
        "value": baby_model.k_wall.magnitude,
        "unit": str(baby_model.k_wall.units),
    },
    "cumulative_tritium_release": {
        label: {
            **{
                form: {
                    "value": gas_stream.get_cumulative_activity(
                        form
                    ).magnitude.tolist(),
                    "unit": str(gas_stream.get_cumulative_activity(form).units),
                }
                for form in ["total", "soluble", "insoluble"]
            },
            "sampling_times": {
                "value": gas_stream.relative_times_as_pint.magnitude.tolist(),
                "unit": str(gas_stream.relative_times_as_pint.units),
            },
        }
        for label, gas_stream in gas_streams.items()
    },
}

# check if the file exists and load it

processed_data_file = "../../data/processed_data.json"

try:
    with open(processed_data_file, "r") as f:
        existing_data = json.load(f)
except FileNotFoundError:
    print(f"Processed data file not found, creating it in {processed_data_file}")
    existing_data = {}

existing_data.update(processed_data)

with open(processed_data_file, "w") as f:
    json.dump(existing_data, f, indent=4)

print(f"Processed data stored in {processed_data_file}")
