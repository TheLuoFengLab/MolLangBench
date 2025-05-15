import argparse
import json
import os
from pathlib import Path
from typing import List, Dict, Any, Optional, Union
import re
from tqdm import tqdm
import pandas as pd
from datasets import load_dataset
from rdkit import Chem
from rdkit.Chem import AllChem

def parse_smiles_from_output(output: str, tag: str = "smiles") -> str:
    """Extract SMILES string from LLM output using specified tags."""
    # Parse the JSON response
    try:
        # Extract the content from the message
        content = output["body"]["choices"][0]["message"]["content"]
    except:
        print("Failed to parse API response, trying to parse raw output")
        content = output
    
    pattern = f"<{tag}>(.*?)</{tag}>"
    match = re.search(pattern, content, re.DOTALL)
    if match:
        return match.group(1).strip()
    return ""

def parse_recognition_output(output: str, result_1_tag: str = "count", result_2_tag: str = "atom_indices") -> Dict[str, Any]:
    """Extract count and atom indices from LLM output using specified tags."""
    # Parse the JSON response
    try:
        # Extract the content from the message
        content = output["body"]["choices"][0]["message"]["content"]
    except:
        print("Failed to parse API response, trying to parse raw output")
        content = output
    
    result = {
        "result_1": None,
        "result_2": None
    }
    
    # Parse result_1
    result_1_pattern = f"<{result_1_tag}>(.*?)</{result_1_tag}>"
    result_1_match = re.search(result_1_pattern, content, re.DOTALL)
    if result_1_match:
        try:
            result["result_1"] = int(result_1_match.group(1).strip())
        except ValueError:
            pass
    
    # Parse result_2
    result_2_pattern = f"<{result_2_tag}>(.*?)</{result_2_tag}>"
    result_2_match = re.search(result_2_pattern, content, re.DOTALL)
    if result_2_match:
        try:
            # Try to parse as list of lists first
            result_2_str = result_2_match.group(1).strip()
            if result_2_str.startswith("[") and result_2_str.endswith("]"):
                # Remove outer brackets and split by inner lists
                inner_lists = result_2_str[1:-1].split("],[")
                result["result_2"] = [
                    [int(idx.strip()) for idx in lst.split(",") if idx.strip()]
                    for lst in inner_lists
                ]
            else:
                # Try parsing as single list
                result["result_2"] = [
                    int(idx.strip()) for idx in result_2_str.split(",") if idx.strip()
                ]
        except (ValueError, SyntaxError):
            pass
    
    return result

def canonicalize_smiles(smiles: str) -> Optional[str]:
    """Canonicalize SMILES string using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None

def evaluate_generation_results(
    results_file: str,
    split: str,
    modality: str
) -> Dict[str, float]:
    """Evaluate generation task results."""
    if modality == "image":
        raise ValueError("Cannot evaluate image outputs for generation task - direct comparison not possible")
        
    # Load dataset
    dataset = load_dataset("ChemFM/MolLangBench", "generation", split=split)
    
    # Load results
    with open(results_file, 'r') as f:
        results = [json.loads(line) for line in f]
    
    # Initialize metrics
    total = len(results)
    valid_smiles = 0
    matching_smiles = 0
    
    # Evaluate each result
    for result in tqdm(results, desc="Evaluating generation results"):
        # Get ground truth
        idx = result.get("custom_id", "").split("_")[-1]
        if not idx.isdigit():
            continue
        idx = int(idx)
        if idx >= len(dataset):
            continue
            
        ground_truth = dataset[idx]["smiles"]
        ground_truth_canon = canonicalize_smiles(ground_truth)
        
        # Get predicted SMILES
        predicted = parse_smiles_from_output(result["response"])
        predicted_canon = canonicalize_smiles(predicted)
        print(predicted_canon, ground_truth_canon)
        
        # Update metrics
        if predicted_canon is not None:
            valid_smiles += 1
            if predicted_canon == ground_truth_canon:
                matching_smiles += 1
    
    return {
        "validity_rate": valid_smiles / total if total > 0 else 0,
        "accuracy": matching_smiles / total if total > 0 else 0
    }

def evaluate_editing_results(
    results_file: str,
    split: str,
    modality: str
) -> Dict[str, float]:
    """Evaluate editing task results."""
    if modality == "image":
        raise ValueError("Cannot evaluate image outputs for editing task - direct comparison not possible")
        
    # Load dataset
    dataset = load_dataset("ChemFM/MolLangBench", "edit", split=split)
    
    # Load results
    with open(results_file, 'r') as f:
        results = [json.loads(line) for line in f]
    
    # Initialize metrics
    total = len(results)
    valid_smiles = 0
    matching_smiles = 0
    
    # Evaluate each result
    for result in tqdm(results, desc="Evaluating editing results"):
        # Get ground truth
        idx = result.get("custom_id", "").split("_")[-1]
        if not idx.isdigit():
            continue
        idx = int(idx)
        if idx >= len(dataset):
            continue
            
        ground_truth = dataset[idx]["edited_smiles"]
        ground_truth_canon = canonicalize_smiles(ground_truth)
        
        # Get predicted SMILES
        predicted = parse_smiles_from_output(result["response"])
        predicted_canon = canonicalize_smiles(predicted)
        
        # Update metrics
        if predicted_canon is not None:
            valid_smiles += 1
            if predicted_canon == ground_truth_canon:
                matching_smiles += 1
    
    return {
        "validity_rate": valid_smiles / total if total > 0 else 0,
        "accuracy": matching_smiles / total if total > 0 else 0
    }

def evaluate_recognition_results(
    results_file: str,
    split: str,
    subtask: str,
    modality: str,
    result_1_tag: str = "count",
    result_2_tag: str = "atom_indices"
) -> Dict[str, float]:
    """Evaluate recognition task results."""
    # Load dataset
    dataset = load_dataset("ChemFM/MolLangBench", f"recognition", split=split)

    # we need to filter the dataset to only include the subtask
    dataset = dataset.filter(lambda x: x["task"] == subtask)
    
    # Load results
    with open(results_file, 'r') as f:
        results = [json.loads(line) for line in f]
    
    # Initialize metrics
    total = len(results)
    result_1_correct = 0
    result_2_correct = 0
    both_correct = 0
    total_result_2 = 0  # Count of non-empty ground truth result_2
    
    # Evaluate each result
    for result in tqdm(results, desc="Evaluating recognition results"):
        # Get ground truth
        idx = result.get("custom_id", "").split("_")[-1]
        if not idx.isdigit():
            continue
        idx = int(idx)
        if idx >= len(dataset):
            continue
            
        # Convert ground truth to appropriate types
        ground_truth_result_1 = int(dataset[idx]["result_1"])
        ground_truth_result_2_str = dataset[idx]["result_2"]
        try:
            ground_truth_result_2 = eval(ground_truth_result_2_str) if ground_truth_result_2_str else []
        except:
            print(f"Warning: Could not parse ground truth result_2: {ground_truth_result_2_str}")
            continue
            
        # Get predicted values
        predicted = parse_recognition_output(
            result["response"],
            result_1_tag=result_1_tag,
            result_2_tag=result_2_tag
        )
        
        # Compare result_1
        result_1_match = predicted["result_1"] == ground_truth_result_1
        
        # Compare result_2
        result_2_match = False
        if predicted["result_2"] is not None and ground_truth_result_2:  # Only compare if ground truth is not empty
            total_result_2 += 1  # Increment counter for non-empty ground truth
            if isinstance(ground_truth_result_2[0], list):
                # Compare list of lists
                if len(predicted["result_2"]) == len(ground_truth_result_2):
                    # we do sort the each individual in the list then sort the list of lists
                    sorted_predicted = sorted([sorted(sub_list) for sub_list in predicted["result_2"]])
                    sorted_ground_truth = sorted([sorted(sub_list) for sub_list in ground_truth_result_2])
                    result_2_match = sorted_predicted == sorted_ground_truth
            else:
                # Compare single list
                result_2_match = sorted(predicted["result_2"]) == sorted(ground_truth_result_2)
        elif not ground_truth_result_2:  # If ground truth is empty, consider result_2 as correct
            result_2_match = True
        
        # Update metrics
        if result_1_match:
            result_1_correct += 1
        if result_2_match:
            result_2_correct += 1
        if result_1_match and result_2_match:
            both_correct += 1
    
    return {
        "result_1_accuracy": result_1_correct / total if total > 0 else 0,
        "result_2_accuracy": result_2_correct / total_result_2 if total_result_2 > 0 else None,
        "both_accuracy": both_correct / total if total > 0 else 0
    }

def main():
    parser = argparse.ArgumentParser(description="Evaluate LLM outputs against ground truth")
    parser.add_argument("--results_file", type=str, required=True,
                      help="Path to JSONL file containing LLM outputs")
    parser.add_argument("--split", type=str, default="test",
                      help="Dataset split to use (default: test)")
    parser.add_argument("--task_type", type=str, required=True, choices=["generation", "editing", "recognition"],
                      help="Task type to evaluate")
    parser.add_argument("--subtask", type=str,
                      help="Subtask to evaluate (only required for recognition task)")
    parser.add_argument("--modality", type=str, required=True,
                      help="Modality to evaluate")
    parser.add_argument("--result_1_tag", type=str, default="count",
                      help="Tag for result_1 in recognition task (default: count)")
    parser.add_argument("--result_2_tag", type=str, default="atom_indices",
                      help="Tag for result_2 in recognition task (default: atom_indices)")
    
    args = parser.parse_args()
    
    # Evaluate based on task type
    if args.task_type == "generation":
        metrics = evaluate_generation_results(
            args.results_file,
            args.split,
            args.modality
        )
        print("\nGeneration Task Metrics:")
        print(f"Validity Rate: {metrics['validity_rate']:.2%}")
        print(f"Accuracy: {metrics['accuracy']:.2%}")
        
    elif args.task_type == "editing":
        metrics = evaluate_editing_results(
            args.results_file,
            args.split,
            args.modality
        )
        print("\nEditing Task Metrics:")
        print(f"Validity Rate: {metrics['validity_rate']:.2%}")
        print(f"Accuracy: {metrics['accuracy']:.2%}")
        
    else:  # recognition
        if not args.subtask:
            raise ValueError("--subtask is required for recognition task")
            
        metrics = evaluate_recognition_results(
            args.results_file,
            args.split,
            args.subtask,
            args.modality,
            args.result_1_tag,
            args.result_2_tag
        )
        print("\nRecognition Task Metrics:")
        print(f"Result 1 Accuracy: {metrics['result_1_accuracy']:.2%}")
        if metrics['result_2_accuracy'] is not None:
            print(f"Result 2 Accuracy: {metrics['result_2_accuracy']:.2%}")
        print(f"Both Correct: {metrics['both_accuracy']:.2%}")

if __name__ == "__main__":
    main() 