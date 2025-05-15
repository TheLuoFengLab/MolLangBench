import argparse
import json
import os
import time
from pathlib import Path
from typing import Optional, Dict, Any
from openai import OpenAI
from tqdm import tqdm

def get_api_key(api_key: Optional[str] = None) -> str:
    """
    Get OpenAI API key from command line argument or environment variable.
    
    Args:
        api_key: API key from command line argument
    
    Returns:
        OpenAI API key
    
    Raises:
        ValueError: If API key is not provided in either command line or environment
    """
    if api_key:
        return api_key
    
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        raise ValueError(
            "OpenAI API key not found. Please provide it either:\n"
            "1. As a command line argument: --api_key your_api_key\n"
            "2. As an environment variable: export OPENAI_API_KEY=your_api_key"
        )
    
    return api_key

def submit_batch_job(
    jobs_file: str,
    api_key: Optional[str] = None,
    organization: Optional[str] = None,
    max_retries: int = 3,
    retry_delay: int = 5
) -> Dict[str, Any]:
    """
    Submit a batch of jobs to OpenAI API.
    
    Args:
        jobs_file: Path to JSONL file containing jobs
        api_key: OpenAI API key (optional if set in environment)
        organization: Optional OpenAI organization ID
        max_retries: Maximum number of retries for failed requests
        retry_delay: Delay between retries in seconds
    
    Returns:
        Dictionary containing batch job information
    """
    # Initialize OpenAI client
    client = OpenAI(api_key=get_api_key(api_key), organization=organization)
    
    for attempt in range(max_retries):
        try:
            # Upload file
            with open(jobs_file, "rb") as f:
                batch_input_file = client.files.create(
                    file=f,
                    purpose="batch"
                )
            
            # Create batch job
            batch_job = client.batches.create(
                input_file_id=batch_input_file.id,
                endpoint="/v1/chat/completions",
                completion_window="24h",
                metadata={
                    "description": f"Batch job for {os.path.basename(jobs_file)}"
                }
            )
            
            return batch_job
            
        except Exception as e:
            if attempt == max_retries - 1:
                raise Exception(f"Failed to submit batch job after {max_retries} attempts: {str(e)}")
            print(f"Attempt {attempt + 1} failed, retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)

def check_batch_status(
    batch_id: str,
    api_key: Optional[str] = None,
    organization: Optional[str] = None,
    max_retries: int = 3,
    retry_delay: int = 5
) -> Dict[str, Any]:
    """
    Check the status of a batch job.
    
    Args:
        batch_id: ID of the batch job
        api_key: OpenAI API key (optional if set in environment)
        organization: Optional OpenAI organization ID
        max_retries: Maximum number of retries for failed requests
        retry_delay: Delay between retries in seconds
    
    Returns:
        Dictionary containing batch job status
    """
    client = OpenAI(api_key=get_api_key(api_key), organization=organization)
    
    for attempt in range(max_retries):
        try:
            batch_job = client.batches.retrieve(batch_id)
            return batch_job
            
        except Exception as e:
            if attempt == max_retries - 1:
                raise Exception(f"Failed to check batch status after {max_retries} attempts: {str(e)}")
            print(f"Attempt {attempt + 1} failed, retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)

def download_batch_results(
    batch_id: str,
    output_file: str,
    api_key: Optional[str] = None,
    organization: Optional[str] = None,
    max_retries: int = 3,
    retry_delay: int = 5
) -> None:
    """
    Download results from a completed batch job.
    
    Args:
        batch_id: ID of the batch job
        output_file: Path to save results
        api_key: OpenAI API key (optional if set in environment)
        organization: Optional OpenAI organization ID
        max_retries: Maximum number of retries for failed requests
        retry_delay: Delay between retries in seconds
    """
    client = OpenAI(api_key=get_api_key(api_key), organization=organization)
    
    for attempt in range(max_retries):
        try:
            # Get batch job to find output file ID
            batch_job = client.batches.retrieve(batch_id)
            output_file_id = batch_job.output_file_id
            
            # Download results
            response = client.files.content(output_file_id)
            
            # Save results
            with open(output_file, "wb") as f:
                f.write(response.content)
            
            return
            
        except Exception as e:
            if attempt == max_retries - 1:
                raise Exception(f"Failed to download batch results after {max_retries} attempts: {str(e)}")
            print(f"Attempt {attempt + 1} failed, retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)

def main():
    parser = argparse.ArgumentParser(description="Submit OpenAI batch jobs and track results")
    subparsers = parser.add_subparsers(dest="mode", help="Operation mode", required=True)
    
    # Submit mode
    submit_parser = subparsers.add_parser("submit", help="Submit a new batch job")
    submit_parser.add_argument("--jobs_file", type=str, required=True,
                            help="Path to JSONL file containing jobs")
    submit_parser.add_argument("--api_key", type=str,
                            help="OpenAI API key (optional if set in OPENAI_API_KEY environment variable)")
    submit_parser.add_argument("--organization", type=str,
                            help="OpenAI organization ID (optional)")
    
    # Retrieve mode
    retrieve_parser = subparsers.add_parser("retrieve", help="Retrieve results from a batch job")
    retrieve_parser.add_argument("--batch_id", type=str, required=True,
                              help="ID of the batch job to retrieve")
    retrieve_parser.add_argument("--output_file", type=str, required=True,
                              help="Path to save batch results")
    retrieve_parser.add_argument("--api_key", type=str,
                              help="OpenAI API key (optional if set in OPENAI_API_KEY environment variable)")
    retrieve_parser.add_argument("--organization", type=str,
                              help="OpenAI organization ID (optional)")
    retrieve_parser.add_argument("--check_interval", type=int, default=60,
                              help="Interval in seconds to check batch status (default: 60)")
    
    args = parser.parse_args()
    
    if args.mode == "submit":
        # Submit batch job
        print("Submitting batch job...")
        batch_info = submit_batch_job(
            jobs_file=args.jobs_file,
            api_key=args.api_key,
            organization=args.organization
        )
        batch_id = batch_info.id
        print(f"\nBatch job submitted successfully!")
        print(f"Batch ID: {batch_id}")
        print(f"Status: {batch_info.status}")
        print("\nTo retrieve results when the job is complete, run:")
        print(f"python scripts/submit_openai_jobs.py retrieve --batch_id {batch_id} --output_file <path_to_save_results>")
        if args.organization:
            print(f"   --organization {args.organization}")
        print("Note: API key can be provided via --api_key or OPENAI_API_KEY environment variable")
            

        
    else:  # retrieve mode
        # Monitor batch status
        print("\nMonitoring batch status...")
        while True:
            status = check_batch_status(
                batch_id=args.batch_id,
                api_key=args.api_key,
                organization=args.organization
            )
            
            print(f"\nStatus: {status.status}")
            if status.status == 'completed':
                print("Batch job completed!")
                break
            elif status.status == 'failed':
                raise Exception(f"Batch job failed: {status.error if hasattr(status, 'error') else 'Unknown error'}")
            
            print(f"Waiting {args.check_interval} seconds before next check...")
            time.sleep(args.check_interval)
        
        # Download results
        print("\nDownloading results...")
        download_batch_results(
            batch_id=args.batch_id,
            output_file=args.output_file,
            api_key=args.api_key,
            organization=args.organization
        )
        print(f"Results saved to: {args.output_file}")

if __name__ == "__main__":
    main() 