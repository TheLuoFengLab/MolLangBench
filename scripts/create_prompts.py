import argparse
import json
import os
import base64
from pathlib import Path
from datasets import load_dataset
from typing import Literal, List, Union, Optional
import ast
from PIL import Image
import io
from tqdm import tqdm

def encode_image(image: Union[str, Image.Image], save_path: Optional[str] = None) -> str:
    """
    Encode image to base64 string and optionally save it.
    
    Args:
        image: Either a path to an image file or a PIL Image object
        save_path: Optional path to save the image
    
    Returns:
        Base64 encoded string of the image
    """
    if isinstance(image, str):
        # If image is a file path
        with open(image, "rb") as image_file:
            image_data = image_file.read()
            if save_path:
                # Copy the file to save_path
                os.makedirs(os.path.dirname(save_path), exist_ok=True)
                with open(save_path, "wb") as f:
                    f.write(image_data)
            return base64.b64encode(image_data).decode('utf-8')
    else:
        # If image is a PIL Image object
        buffered = io.BytesIO()
        image.save(buffered, format="PNG")
        image_data = buffered.getvalue()
        if save_path:
            # Save the image to save_path
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            with open(save_path, "wb") as f:
                f.write(image_data)
        return base64.b64encode(image_data).decode('utf-8')

def parse_target_atoms(target_atoms_str: str) -> Union[int, List[int], None]:
    """
    Parse target atoms string into either a single integer, a list of two integers, or None.
    Returns None for empty or null strings.
    """
    if not target_atoms_str or target_atoms_str.strip() == "":
        return None
        
    try:
        # Try to parse as a list first
        atoms = ast.literal_eval(target_atoms_str)
        if isinstance(atoms, list) and len(atoms) == 2:
            return atoms
        elif isinstance(atoms, int):
            return atoms
        else:
            raise ValueError(f"Invalid target atoms format: {target_atoms_str}")
    except (ValueError, SyntaxError):
        # If not a list, try to parse as a single integer
        try:
            return int(target_atoms_str)
        except ValueError:
            raise ValueError(f"Invalid target atoms format: {target_atoms_str}")

def create_generation_prompts(
    output_file: str,
    modality: Literal["image", "smiles"],
    split: str = "test",
    prompt_template_dir: str = "prompts/generation"
) -> None:
    """
    Create JSONL file containing prompts for the generation task.
    
    Args:
        output_file: Path to output JSONL file
        modality: Either "image" or "smiles"
        split: Dataset split to use (default: "test")
        prompt_template_dir: Directory containing prompt templates
    """
    if split != "test":
        raise ValueError("Only test split is available for generation task")
    
    # Load dataset with generation config
    dataset = load_dataset("ChemFM/MolLangBench", "generation", split=split)
    
    # Load prompt template
    template_path = Path(prompt_template_dir) / f"prompt_{modality}.txt"
    if not template_path.exists():
        raise FileNotFoundError(f"Prompt template not found at {template_path}")
    
    with open(template_path, "r") as f:
        template = f.read().strip()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Generate prompts and write to JSONL
    with open(output_file, "w") as f:
        for idx, item in tqdm(enumerate(dataset), total=len(dataset), desc="Creating generation prompts"):
            # Format prompt using template
            prompt = template.format(structure_description=item["structure_description"])
            id = f"generation_{split}_{modality}_{idx}"
            
            # Create output entry
            output_entry = {
                "id": id,
                "prompt": prompt
            }
            
            # Write to JSONL
            f.write(json.dumps(output_entry) + "\n")

def create_editing_prompts(
    output_file: str,
    modality: Literal["image", "smiles"],
    split: str = "test",
    prompt_template_dir: str = "prompts/editing",
    save_images: bool = False
) -> None:
    """
    Create JSONL file containing prompts for the editing task.
    
    Args:
        output_file: Path to output JSONL file
        modality: Either "image" or "smiles"
        split: Dataset split to use (default: "test")
        prompt_template_dir: Directory containing prompt templates
        save_images: Whether to save input images (only for image modality)
    """
    if split != "test":
        raise ValueError("Only test split is available for editing task")
    
    # Load dataset with editing config
    dataset = load_dataset("ChemFM/MolLangBench", "edit", split=split)
    
    # Load prompt template
    template_path = Path(prompt_template_dir) / f"prompt_{modality}.txt"
    if not template_path.exists():
        raise FileNotFoundError(f"Prompt template not found at {template_path}")
    
    with open(template_path, "r") as f:
        template = f.read().strip()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Generate prompts and write to JSONL
    with open(output_file, "w") as f:
        for idx, item in tqdm(enumerate(dataset), total=len(dataset), desc="Creating editing prompts"):
            # Format prompt using template
            if modality == "smiles":
                prompt = template.format(
                    original_smiles=item["original_smiles"],
                    edit_instructions=item["edit_instructions"]
                )
                output_entry = {
                    "id": f"editing_{split}_{modality}_{idx}",
                    "prompt": prompt
                }
            else:  # image modality
                prompt = template.format(edit_instructions=item["edit_instructions"])
                # Encode the image
                image_save_path = os.path.join(os.path.dirname(output_file), "images", f"editing_{split}_{modality}_{idx}.png") if save_images else None
                encoded_image = encode_image(item["original_image"], image_save_path)
                output_entry = {
                    "id": f"editing_{split}_{modality}_{idx}",
                    "prompt": prompt,
                    "input_image": encoded_image
                }
            
            # Write to JSONL
            f.write(json.dumps(output_entry) + "\n")

def create_recognition_prompts(
    output_file: str,
    modality: Literal["image", "smiles"],
    subtask: str,
    split: str = "test",
    prompt_template_dir: str = "prompts/recognition",
    save_images: bool = False
) -> None:
    """
    Create JSONL file containing prompts for the recognition task.
    
    Args:
        output_file: Path to output JSONL file
        modality: Either "image" or "smiles"
        subtask: The specific recognition subtask
        split: Dataset split to use (default: "test")
        prompt_template_dir: Directory containing prompt templates
        save_images: Whether to save input images (only for image modality)
    """
    # Load dataset with recognition config and subtask
    dataset = load_dataset("ChemFM/MolLangBench", f"recognition", split=split)
    
    # we need to filter the dataset to only include the subtask
    dataset = dataset.filter(lambda x: x["task"] == subtask)
    
    # Load prompt template
    template_path = Path(prompt_template_dir) / subtask / f"prompt_{modality}.txt"
    if not template_path.exists():
        raise FileNotFoundError(f"Prompt template not found at {template_path}")
    
    with open(template_path, "r") as f:
        template = f.read().strip()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Generate prompts and write to JSONL
    with open(output_file, "w") as f:
        for idx, item in tqdm(enumerate(dataset), total=len(dataset), desc=f"Creating recognition prompts for {subtask}"):
            # Format prompt using template
            if modality == "smiles":
                # Handle target atoms if present
                if "target_atoms" in item:
                    target_atoms = parse_target_atoms(item["target_atoms"])
                    if target_atoms is None:
                        # No target atoms specified
                        prompt = template.format(smiles=item["smiles"])
                    elif isinstance(target_atoms, list):
                        prompt = template.format(
                            smiles=item["smiles"],
                            target_atom_0=target_atoms[0],
                            target_atom_1=target_atoms[1]
                        )
                    else:
                        prompt = template.format(
                            smiles=item["smiles"],
                            target_atom=target_atoms
                        )
                else:
                    prompt = template.format(smiles=item["smiles"])
                
                output_entry = {
                    "id": f"recognition_{subtask}_{split}_{modality}_{idx}",
                    "prompt": prompt
                }
            else:  # image modality
                prompt = template.format()
                # Encode the image
                image_save_path = os.path.join(os.path.dirname(output_file), "images", f"recognition_{subtask}_{split}_{modality}_{idx}.png") if save_images else None
                encoded_image = encode_image(item["image"], image_save_path)
                output_entry = {
                    "id": f"recognition_{subtask}_{split}_{modality}_{idx}",
                    "prompt": prompt,
                    "input_image": encoded_image
                }
            
            # Write to JSONL
            f.write(json.dumps(output_entry) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Create JSONL files for prompts")
    parser.add_argument("--task_type", type=str, required=True, choices=["recognition", "generation", "editing"],
                      help="Type of task (recognition, generation, or editing)")
    parser.add_argument("--modality", type=str, required=True, choices=["image", "smiles"],
                      help="Modality for molecule representation (image or smiles)")
    parser.add_argument("--output_file", type=str, required=True,
                      help="Path to output JSONL file")
    parser.add_argument("--split", type=str, default="test",
                      help="Dataset split to use (default: test)")
    parser.add_argument("--recognition_subtask", type=str,
                      help="Subtask type for recognition task (required only for recognition task)")
    parser.add_argument("--save_images", action="store_true",
                      help="Whether to save input images (only for image modality in recognition and editing tasks)")
    
    args = parser.parse_args()
    
    if args.task_type == "generation":
        create_generation_prompts(
            output_file=args.output_file,
            modality=args.modality,
            split=args.split
        )
    elif args.task_type == "recognition":
        if not args.recognition_subtask:
            raise ValueError("recognition_subtask is required for recognition task")
        create_recognition_prompts(
            output_file=args.output_file,
            modality=args.modality,
            subtask=args.recognition_subtask,
            split=args.split,
            save_images=args.save_images
        )
    elif args.task_type == "editing":
        create_editing_prompts(
            output_file=args.output_file,
            modality=args.modality,
            split=args.split,
            save_images=args.save_images
        )

if __name__ == "__main__":
    main() 