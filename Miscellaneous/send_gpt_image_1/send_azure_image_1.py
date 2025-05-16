import os  
import base64
from openai import AzureOpenAI  
import argparse
from tqdm import tqdm
import json
from mimetypes import guess_type
from PIL import Image
from io import BytesIO
import requests

arg_parser = argparse.ArgumentParser()
# define the start and end index
arg_parser.add_argument("--start", type=int, default=0, help="The start index of the requests")
arg_parser.add_argument("--end", type=int, default=100, help="The end index of the requests")
args = arg_parser.parse_args()

endpoint = os.getenv("AZURE_OPENAI_ENDPOINT", "AZURE_OPENAI_ENDPOINT")  # Replace with your actual endpoint
deployment = os.getenv("DEPLOYMENT_NAME", "gpt-image-1")
api_version = os.getenv("OPENAI_API_VERSION", "2025-04-01-preview")
subscription_key = os.getenv("AZURE_OPENAI_API_KEY", "API_KEY")  # Replace with your actual API key

def decode_and_save_image(b64_data, output_filename):
  image = Image.open(BytesIO(base64.b64decode(b64_data)))
  #image.show()
  image.save(output_filename)

def save_all_images_from_response(response_data, filename_prefix):
  for idx, item in enumerate(response_data['data']):
    b64_img = item['b64_json']
    filename = f"{filename_prefix}_{idx+1}.png"
    decode_and_save_image(b64_img, filename)

base_path = f'openai/deployments/{deployment}/images'
params = f'?api-version={api_version}'

# read the jsonl file
with open("batch_request.jsonl", "r") as f:
    prompts = f.readlines()

prompts = prompts[args.start:args.end]

# open the output jsonl file
json_file_name = f"output_{args.start}_{args.end}.jsonl"

existing_custom_ids = set()
if os.path.exists(json_file_name):
    with open(json_file_name, "r") as f:
        for line in f:
            try:
                entry = json.loads(line)
                existing_custom_ids.add(entry["custom_id"])
            except json.JSONDecodeError:
                continue  # Skip any corrupted lines

output_dir = f"image_output"
os.makedirs(output_dir, exist_ok=True)

json_file = open(json_file_name, "a")
for request in tqdm(prompts):
    request_dict = json.loads(request)
    custom_id = request_dict.get("custom_id")
    if custom_id in existing_custom_ids:
        print(f"Skipping custom_id {custom_id}, already processed.")
        continue

    prompt = request_dict['body']['messages'][0]['content'][0]['text']
    image_path = request_dict['image_path']

    # Check if output image already exists
    output_filename = f"{output_dir}/{custom_id}_1.png"
    if os.path.exists(output_filename):
        print(f"Skipping custom_id {custom_id}, output image already exists.")
        continue

    edit_url = f"{endpoint}{base_path}/edits{params}"
    edit_body = {
      "prompt": prompt,
      "n": 1,
      "size": "1024x1024",
      "quality": "medium"
    }

    files = {
      "image": (image_path, open(image_path, "rb"), "image/png"),
    }
    output_filename = f"{output_dir}/{custom_id}"
    try:
        edit_response = requests.post(
          edit_url,
          headers={'Api-Key': subscription_key},
          data=edit_body,
          files=files
        ).json()
        save_all_images_from_response(edit_response, output_filename)
    except Exception as e:
        print(f"An error occurred: {e}")
