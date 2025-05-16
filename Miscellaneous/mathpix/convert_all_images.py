import os
import json
import requests
import time
from pathlib import Path
from tqdm import tqdm

def get_processed_images(output_file):
    """Get set of already processed image paths from output file."""
    processed = set()
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            for line in f:
                try:
                    data = json.loads(line.strip())
                    if 'image_path' in data:
                        processed.add(data['image_path'])
                except json.JSONDecodeError:
                    continue
    return processed

def convert_image_to_smiles(image_path, app_id, app_key):
    """Convert a single image to SMILES using Mathpix API."""
    try:
        response = requests.post(
            "https://api.mathpix.com/v3/text",
            files={"file": open(image_path, "rb")},
            data={
                "options_json": json.dumps({
                    "rm_spaces": True,
                    "include_line_data": True,
                    "include_smiles": True,
                })
            },
            headers={
                "app_id": app_id,
                "app_key": app_key,
            }
        )
        result = response.json()
        # add the image_path to the result
        result["image_path"] = str(image_path)
        return result
    except Exception as e:
        return {
            "image_path": str(image_path),
            "smiles": "",
            "error": str(e)
        }

def process_directory(directory_path, output_file, app_id, app_key):
    """Process all PNG images in a directory and save results to JSONL."""
    directory = Path(directory_path)
    image_files = sorted([f for f in directory.glob("*.png")])
    
    # Get set of already processed images
    processed_images = get_processed_images(output_file)
    
    # Create a dictionary to store results in order
    results_dict = {}
    
    # First, load existing results to maintain order
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            for line in f:
                try:
                    data = json.loads(line.strip())
                    if 'image_path' in data:
                        results_dict[data['image_path']] = data
                except json.JSONDecodeError:
                    continue
    
    # Process new images
    new_images = [f for f in image_files if str(f) not in processed_images]
    if not new_images:
        print(f"No new images to process in {directory.name}")
        return
    
    for image_file in tqdm(new_images, desc=f"Processing {directory.name}"):
        result = convert_image_to_smiles(image_file, app_id, app_key)
        results_dict[str(image_file)] = result
        
        # Add 10-second delay between requests
        if image_file != new_images[-1]:  # Don't delay after the last image
            print(f"\nWaiting 10 seconds before next request...")
            time.sleep(10)
    
    # Write all results in order
    with open(output_file, 'w') as f:
        for image_file in image_files:
            if str(image_file) in results_dict:
                f.write(json.dumps(results_dict[str(image_file)]) + "\n")

def main():
    # Mathpix API credentials
    APP_ID = "APP_ID"  # Replace with your actual APP_ID
    APP_KEY = "API_KEY" # Replace with your actual APP_KEY
    
    # Define directories and output files
    gpt_image_dir = Path("../gpt_image_1/image_output")
    
    output_gpt = "gpt_image_smiles.jsonl"
    
    # Process both directories
    print("Processing GPT image directory...")
    process_directory(gpt_image_dir, output_gpt, APP_ID, APP_KEY)
    
    print("\nConversion complete! Results saved to:")
    print(f"- {output_gpt}")

if __name__ == "__main__":
    main() 