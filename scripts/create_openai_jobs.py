import argparse
import json
from pathlib import Path
from typing import Optional

def create_openai_jobs(
    prompt_file: str,
    output_file: str,
    model: str,
    custom_id_prefix: Optional[str] = None
) -> None:
    """
    Create JSONL file for OpenAI API submission.
    
    Args:
        prompt_file: Path to the input prompt JSONL file
        output_file: Path to output JSONL file for OpenAI API
        model: OpenAI model to use (e.g. o1-mini, o1, o3)
        custom_id_prefix: Optional prefix for custom IDs
    """
    # Read prompts from input file
    with open(prompt_file, "r") as f:
        prompts = [json.loads(line) for line in f]
    
    # Create output directory if it doesn't exist
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    # Generate OpenAI API jobs and write to JSONL
    with open(output_file, "w") as f:
        for prompt in prompts:
            # Create custom ID if prefix is provided
            custom_id = f"{custom_id_prefix}_{prompt['id']}" if custom_id_prefix else prompt['id']
            
            # Check if prompt has an input_image field
            if "input_image" in prompt:
                # Create message with both text and image
                messages = [
                    {
                        "type": "text",
                        "text": prompt["prompt"]
                    },
                    {
                        "type": "image_url",
                        "image_url": {
                            "url": f"data:image/png;base64,{prompt['input_image']}"
                        }
                    }
                ]
            else:
                # Regular text prompt
                messages = [{
                    "type": "text",
                    "text": prompt["prompt"]
                }]
            
            # Create OpenAI API job entry
            job_entry = {
                "custom_id": custom_id,
                "method": "POST",
                "url": "/v1/chat/completions",
                "body": {
                    "model": model,
                    "messages": [
                        {
                            "role": "user",
                            "content": messages
                        }
                    ]
                }
            }
            
            # Write to JSONL
            f.write(json.dumps(job_entry) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Create JSONL files for OpenAI API submission")
    parser.add_argument("--prompt_file", type=str, required=True,
                      help="Path to input prompt JSONL file")
    parser.add_argument("--output_file", type=str, required=True,
                      help="Path to output JSONL file for OpenAI API")
    parser.add_argument("--model", type=str, required=True,
                      help="OpenAI model to use (e.g. o1-mini, o1, o3)")
    parser.add_argument("--custom_id_prefix", type=str, default=None,
                      help="Optional prefix for custom IDs")
    
    args = parser.parse_args()
    
    create_openai_jobs(
        prompt_file=args.prompt_file,
        output_file=args.output_file,
        model=args.model,
        custom_id_prefix=args.custom_id_prefix
    )

if __name__ == "__main__":
    main() 