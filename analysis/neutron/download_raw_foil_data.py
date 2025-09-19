from pathlib import Path
import zipfile
import requests


def download_and_extract_foil_data(url: str, extracted_path: Path):

    output_filepath = Path("../../data/neutron_detection/foil_data.zip")

    if extracted_path.exists():
        print(f"Directory already exists: {extracted_path}")
    else:
        # URL of the file

        # Download the file
        print(f"Downloading data from {url}...")
        response = requests.get(url)
        if response.status_code == 200:
            print("Download successful!")
            # Save the file to the specified directory
            with open(output_filepath, "wb") as f:
                f.write(response.content)
            print(f"File saved to: {output_filepath}")
        else:
            print(f"Failed to download file. HTTP Status Code: {response.status_code}")

        # Extract the zip file

        # Ensure the extraction directory exists
        extracted_path.mkdir(parents=True, exist_ok=True)

        # Unzip the file
        with zipfile.ZipFile(output_filepath, "r") as zip_ref:
            zip_ref.extractall(extracted_path)
        print(f"Files extracted to: {extracted_path}")

        # Delete the zip file after extraction
        output_filepath.unlink(missing_ok=True)

