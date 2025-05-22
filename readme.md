<a id="readme-top"></a>

<!-- PROJECT SHIELDS -->
[stars-shield]: https://img.shields.io/github/stars/TheLuoFengLab/MolLangBench.svg?style=flat-square&color=b75347
[stars-url]: https://github.com/TheLuoFengLab/MolLangBench/stargazers
[forks-shield]: https://img.shields.io/github/forks/TheLuoFengLab/MolLangBench.svg?style=flat-square&color=df7e66
[forks-url]: https://github.com/TheLuoFengLab/MolLangBench/network/members
[issues-shield]: https://img.shields.io/github/issues/TheLuoFengLab/MolLangBench.svg?style=flat-square&color=edc775
[issues-url]: https://github.com/TheLuoFengLab/MolLangBench/issues
[license-shield]: https://img.shields.io/badge/License-MIT-lightgrey.svg?style=flat-square&color=94b594
[license-url]: https://github.com/TheLuoFengLab/MolLangBench/blob/main/LICENSE

<!-- PROJECT LOGO -->
<br />
<div align="center">
  <h1 align="center">MolLangBench: A Comprehensive Benchmark for Language-Prompted Molecular Structure Recognition, Editing, and Generation</h1>

  [![Stargazers][stars-shield]][stars-url]
  [![GitHub License][license-shield]][license-url]

  <p align="center">
    <a href="https://arxiv.org/abs/2505.15054">
      <img src="https://info.arxiv.org/brand/images/brand-supergraphic.jpg" alt="arxiv" width="25" height="25" style="vertical-align: middle; margin-right: 0px;">
    </a>    
    <a href="https://arxiv.org/abs/2505.15054">
      ArXiv
    </a>
    |
    <a href="https://huggingface.co/datasets/ChemFM/MolLangBench">
      <img src="https://huggingface.co/front/assets/huggingface_logo.svg" alt="Hugging Face" width="20" height="20">
    </a>
    <a href="https://huggingface.co/datasets/ChemFM/MolLangBench">Hugging Face Dataset</a>
    |
    <a href="https://discord.gg/hpW7sdMQGP">
      <img src="https://camo.githubusercontent.com/ae76bfbcd3ea4af324682842213b28d9a7ebdd8791d8531d1b7e3b8b4d2a0302/68747470733a2f2f6564656e742e6769746875622e696f2f537570657254696e7949636f6e732f696d616765732f7376672f646973636f72642e737667" alt="Discord" width="25" height="25">
    </a>
    <a href="https://discord.gg/hpW7sdMQGP">Discord</a>
  </p>

</div>

<details>
  <summary>Table of Contents</summary>
    <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
    </li>
    <li>
      <a href="#quick-start">Quick Start</a>
    </li>
    <li>
      <a href="#usage">Usage</a>
    </li>
    <li><a href="#miscellaneous">Miscellaneous</a></li>
    <li>
      <a href="#benchmark-results">Benchmark Results</a>
    </li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#citation">Citation</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>
<!--
- [Citation](#citation)
-->

## About The Project

**MolLangBench** is a comprehensive benchmark designed to evaluate the fundamental capabilities of AI models in language-prompted molecular structure recognition, editing, and generation.

This repository provides:
- Code and examples to load and use the dataset directly from the [Hugging Face Dataset](https://huggingface.co/datasets/ChemFM/MolLangBench)
- Evaluation scripts and prompt templates to test OpenAI models (e.g., o1, o3, o4-mini) using either molecular images or SMILES strings as inputs

It is straightforward to extend this repository to evaluate other language or multimodal models by adapting the provided input formatting and evaluation templates.

## Getting Started

You can easily set up the required environment using Conda by following these steps:

* Clone the repository

    ```bash
    git clone https://github.com/TheLuoFengLab/MolLangBench.git
    cd MolLangBench
    ```

* Install dependencies
    ```bash
    pip install -r requirements.txt
    ```
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->

## Quick Start


Load the dataset in Python:

```python
from datasets import load_dataset

# Recognition (train + test)
rec_train = load_dataset("ChemFM/MolLangBench", name="recognition", split="train")
rec_test  = load_dataset("ChemFM/MolLangBench", name="recognition", split="test")

# Filter one specific subtask
subtask = "one_hop_neighbors"
subset  = rec_test.filter(lambda x: x["task"] == subtask)

# Editing (test only)
edit = load_dataset("ChemFM/MolLangBench", name="edit", split="test")

# Generation (test only)
gen  = load_dataset("ChemFM/MolLangBench", name="generation", split="test")
```
<p align="right">(<a href="#readme-top">back to top</a>)</p>


<a id="usage-top"></a>


## Usage

We provide end-to-end scripts for:

1. Generating prompt files (`.jsonl`)
2. Creating batch job inputs for the OpenAI API
3. Submitting jobs & retrieving outputs
4. Computing evaluation metrics

Below is a step-by-step example workflow for the **“one hop neighbors”** recognition subtask.

---

### 1. Prepare Prompts

Prompt templates for all tasks and modalities (SMILES and image) are located in the [`prompts`](./prompts/) folder. To generate a `.jsonl` prompt file, run:

```bash
python scripts/create_prompts.py \
    --task_type <recognition|editing|generation> \
    --recognition_subtask <recognition_subtask_name_if_applicable> \
    --modality <smiles|image> \
    --split <train|test> \
    --output_file <output_jsonl_path>
````

**Example:**
For the "one hop neighbors" recognition subtask (SMILES modality, test split):

```bash
python scripts/create_prompts.py \
    --task_type recognition \
    --recognition_subtask one_hop_neighbors \
    --modality smiles \
    --split test \
    --output_file exps/one_hop_neighbors/prompts.jsonl
```

> For image modality, the image is included as a Base64-encoded string in the `.jsonl` file.

---

### 2. Create OpenAI Batch Job File

Generate a batch job file for your prompts and desired model:

```bash
python scripts/create_openai_jobs.py \
    --prompt_file <prompts_jsonl_path> \
    --output_file <batch_job_jsonl_path> \
    --model <model_id> \
    --custom_id_prefix <optional_prefix>
```

**Example:**
Using the `o4-mini` model for "one hop neighbors":

```bash
python scripts/create_openai_jobs.py \
    --prompt_file exps/one_hop_neighbors/prompts.jsonl \
    --output_file exps/one_hop_neighbors/o4-mini/batch_input.jsonl \
    --model o4-mini \
    --custom_id_prefix o4_mini
```

---

### 3. Submit Jobs & Retrieve Outputs

Submit the batch job to the OpenAI API:

```bash
python scripts/submit_openai_jobs.py submit \
    --jobs_file <batch_job_jsonl_path> \
    [--api_key YOUR_API_KEY] \
    [--organization YOUR_ORG_ID]
```

> You can also set your API key as an environment variable. The organization ID is optional.

**Example:**

```bash
python scripts/submit_openai_jobs.py submit \
    --jobs_file exps/one_hop_neighbors/o4-mini/batch_input.jsonl
```

This command will print a `batch_id` for your job.

To retrieve the results (periodically checks until the job is complete):

```bash
python scripts/submit_openai_jobs.py retrieve \
    --batch_id <BATCH_ID> \
    --output_file <results_jsonl_path> \
    [--api_key YOUR_API_KEY] \
    [--organization YOUR_ORG_ID] \
    [--check_interval 60]
```

**Example:**

```bash
python scripts/submit_openai_jobs.py retrieve \
    --batch_id BATCH_ID \
    --output_file exps/one_hop_neighbors/o4-mini/results.jsonl
```

---

### 4. Evaluate the Results

Evaluate the model outputs with:

```bash
python scripts/evaluate_results.py \
    --results_file <results_jsonl_path> \
    --task_type <recognition|editing|generation> \
    --subtask <subtask_name> \
    --modality <smiles|image>
```

**Example:**

```bash
python scripts/evaluate_results.py \
    --results_file exps/one_hop_neighbors/o4-mini/results.jsonl \
    --task_type recognition \
    --subtask one_hop_neighbors \
    --modality smiles
```
This will print out the evaluation metrics for your selected task and model.
> The default result tags are `<count>` and `<atom_indices>`. For certain tasks, you may need to specify custom result tags using the `--result_1_tag <result_1_tag>` argument.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

> **⚠️ Warning:**  
> For **molecule editing and generation tasks using the `gpt-image-1` model**, the OpenAI API does **not** support batch job submissions.  
> - You must submit each image generation or editing request individually using the [OpenAI Images API](https://platform.openai.com/docs/api-reference/images).
> - For automatic evaluation, you will also need a [Mathpix](https://mathpix.com/) account and API key to convert the generated molecular images back to SMILES strings.
>
> Example scripts for these tasks are provided in the [`Miscellaneous`](./Miscellaneous) folder.  
> We plan to provide a more streamlined workflow in the future.

---

## Miscellaneous

The [`Miscellaneous`](./Miscellaneous) folder contains helpful scripts and utilities, including:

1. **Ground Truth Collection**  
   Scripts for collecting ground truth information for each recognition task using RDKit.

2. **Image-to-SMILES Conversion**  
   Scripts to call the Mathpix API for converting molecule images to SMILES strings for automated evaluation.

3. **Per-Image OpenAI API Submission**  
   Scripts to submit image generation and editing requests (for the image modality) to the OpenAI `gpt-image-1` API one-by-one, as batch jobs are not currently supported.

More utilities and improvements will be added in the future.
<p align="right">(<a href="#readme-top">back to top</a>)</p>



## Benchmark Results
### Molecular Structure Recognition

Below are the complete evaluation results for all molecular structure recognition subtasks across a wide range of large language models and vision-language models.  

<details>
<summary><b>Complete Results for Molecular Structure Recognition Tasks</b> (click to expand)</summary>

| Task                | GPT‑4o      | GPT‑4.5-preview | GPT‑4.1     | o1‑mini     | o1          | o3-mini     | DeepSeek‑R1 | R1‑70B      | o3             | o4‑mini     | o3 (image)  | o4-mini (image) |
| ------------------- | ----------- | --------------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ---------------| ----------- | ----------- | --------------- |
| One‑hop neighbors   | 0.355/0.140 | 0.600/0.425     | 0.570/0.330 | 0.735/0.640 | 0.825/0.720 | 0.870/0.820 | 0.825/0.710 | 0.585/0.430 | **0.935/0.895**| 0.880/0.845 | 0.890/0.855 | 0.840/0.780     |
| Two‑hop neighbors   | 0.215/0.055 | 0.280/0.100     | 0.400/0.210 | 0.465/0.350 | 0.745/0.560 | 0.820/0.740 | 0.610/0.475 | 0.305/0.135 | **0.935/0.825**| 0.870/0.790 | 0.770/0.705 | 0.775/0.690     |
| Three‑hop neighbors | 0.165/0.015 | 0.355/0.165     | 0.265/0.140 | 0.400/0.265 | 0.560/0.400 | 0.825/0.705 | 0.550/0.385 | 0.300/0.130 | **0.925/0.830**| 0.775/0.710 | 0.695/0.600 | 0.660/0.575     |
| Quaternary carbons  | 0.530/0.290 | 0.690/0.435     | 0.740/0.440 | 0.615/0.470 | 0.865/0.665 | 0.835/0.740 | 0.780/0.680 | 0.440/0.330 | **0.935/0.865**| 0.845/0.750 | 0.670/0.600 | 0.720/0.665     |
| Ring junctions      | 0.285/0.080 | 0.495/0.185     | 0.485/0.210 | 0.325/0.175 | 0.575/0.470 | 0.580/0.520 | 0.535/0.420 | 0.255/0.160 | **0.685/0.650**| 0.590/0.570 | 0.660/0.595 | 0.615/0.555     |
| Bond connection     | 0.448       | 0.472           | 0.336       | 0.698       | 0.758       | 0.832       | 0.802       | 0.564       | **0.950**      | 0.880       | 0.626       | 0.706           |
| Halogen atoms       | 0.845/0.290 | 0.905/0.420     | 0.900/0.355 | 0.920/0.570 | 0.975/0.740 | 0.955/0.710 | 0.970/0.735 | 0.740/0.375 | **0.965/0.860**| 0.965/0.820 | 0.855/0.815 | 0.920/0.860     |
| Aldehyde            | 0.855/0.570 | 0.965/0.610     | 0.945/0.730 | 0.855/0.725 | 0.970/0.825 | 0.985/0.920 | 0.960/0.835 | 0.715/0.585 | **0.990/0.960**| 0.985/0.945 | 0.925/0.925 | 0.975/0.965     |
| Amide               | 0.505/0.180 | 0.570/0.205     | 0.635/0.315 | 0.585/0.340 | 0.715/0.440 | 0.685/0.510 | 0.635/0.415 | 0.495/0.205 | **0.765/0.650**| 0.755/0.610 | 0.565/0.500 | 0.735/0.665     |
| Carboxyl            | 0.760/0.260 | 0.885/0.235     | 0.900/0.485 | 0.840/0.580 | 0.965/0.675 | 0.955/0.760 | 0.900/0.660 | 0.820/0.495 | **0.985/0.845**| 0.950/0.725 | 0.785/0.750 | 0.870/0.820     |
| Ester               | 0.600/0.145 | 0.760/0.285     | 0.780/0.330 | 0.675/0.325 | 0.935/0.500 | 0.895/0.645 | 0.680/0.400 | 0.615/0.270 | **0.955/0.780**| 0.950/0.640 | 0.720/0.505 | 0.840/0.595     |
| Ketone              | 0.530/0.155 | 0.750/0.260     | 0.870/0.435 | 0.750/0.465 | 0.925/0.600 | 0.985/0.745 | 0.880/0.600 | 0.770/0.370 | **0.985/0.865**| 0.985/0.795 | 0.765/0.675 | 0.850/0.775     |
| Benzene             | 0.490/0.145 | 0.540/0.105     | 0.660/0.155 | 0.530/0.235 | 0.720/0.360 | 0.725/0.565 | 0.595/0.385 | 0.500/0.190 | **0.880/0.695**| 0.730/0.550 | 0.675/0.405 | 0.680/0.485     |
| Furan               | 0.295/0.265 | 0.820/0.325     | 0.905/0.515 | 0.780/0.500 | 0.920/0.660 | 0.865/0.745 | 0.895/0.710 | 0.850/0.490 | **0.975/0.845**| 0.940/0.790 | 0.890/0.820 | 0.870/0.815     |
| Pyridine            | 0.555/0.225 | 0.525/0.250     | 0.730/0.365 | 0.685/0.375 | 0.765/0.555 | 0.860/0.740 | 0.685/0.520 | 0.630/0.340 | **0.925/0.825**| 0.835/0.750 | 0.715/0.585 | 0.790/0.665     |
| Thiophene           | 0.860/0.385 | 0.840/0.325     | 0.880/0.480 | 0.840/0.605 | 0.915/0.690 | 0.940/0.795 | 0.920/0.705 | 0.850/0.565 | **0.970/0.890**| 0.925/0.820 | 0.960/0.855 | 0.920/0.855     |
| Bond stereo         | 0.390       | 0.395           | **0.670**   | 0.425       | 0.330       | 0.310       | 0.310       | 0.345       |   0.480        | 0.325       | 0.575       | 0.640           |
| Chiral stereo       | 0.440       | 0.395           | 0.530       | 0.465       | 0.510       | 0.435       | 0.440       | 0.495       | **0.545**      | 0.520       | 0.510       | 0.495           |
| **Average**         | 0.507/0.249 | 0.625/0.311     | 0.678/0.391 | 0.644/0.456 | 0.776/0.581 | 0.798/0.680 | 0.721/0.566 | 0.571/0.360 | **0.877/0.792**| 0.817/0.713 | 0.736/0.661 | 0.772/0.700     |

---

* Each entry reports recognition accuracy / localization accuracy where applicable.
* Tasks with only recognition evaluation show a single recognition accuracy value.
* **Bold** values indicate the best performance among all evaluated language models.
* "o3 (image)" and "o4-mini (image)" indicate vision-language models evaluated on molecular images.
</details>

---

### Molecule Editing and Generation Benchmark Results

Below are the complete evaluation results for molecule editing and generation tasks across all evaluated language and vision-language models.  

<details>
<summary><b>Complete Results for Molecular Structure Recognition Tasks</b> (click to expand)</summary>

| Task                | GPT‑4o     | GPT‑4.5-preview | GPT‑4.1    | o1‑mini   | o1        | o3-mini   | DeepSeek‑R1 | R1‑70B    | o3        | o4‑mini    | GPT Image 1 |
|---------------------|------------|-----------------|------------|-----------|-----------|-----------|-------------|-----------|-----------|------------|-------------|
| Molecule editing    | 0.725/0.400| 0.950/0.570     | 0.835/0.465| 0.710/0.385| 0.845/0.635| 0.805/0.650| 0.720/0.485 | 0.675/0.375| **0.945/0.785**| 0.920/0.690| 0.135        |
| Molecule generation | 0.525/0.005| 0.800/0.055     | 0.710/0.035| 0.335/0.035| 0.385/0.100| 0.450/0.175| 0.400/0.045 | 0.205/0.010| **0.670/0.290**| 0.600/0.260| 0.000        |

- Each entry reports **SMILES validity / accuracy**.
- **Bold** entries highlight the best performance.
</details>

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- CONTACT -->
## Contact

Main Developer: [Feiyang Cai](mailto:feiyang@clemson.edu) - feiyang@clemson.edu  
Project Supervisor: [Feng Luo](mailto:luofeng@clemson.edu) - luofeng@clemson.edu  

Join our community on [Discord](https://discord.gg/hpW7sdMQGP) to stay updated or ask questions.

<p align="right">(<a href="#readme-top">back to top</a>)</p>


## Citation
If you find our work valuable, please consider giving the project a star and citing it in your research:
```
@article{MolLangBench,
      title={MolLangBench: A Comprehensive Benchmark for Language-Prompted Molecular Structure Recognition, Editing, and Generation},  
      author={Feiyang Cai and Jiahui Bai and Tao Tang and Joshua Luo and Tianyu Zhu and Ling Liu and Feng Luo},
      year={2025},
      journal = {arXiv preprint arXiv:2505.15054},
}
```
Thank you for your support!
<p align="right">(<a href="#top">back to top</a>)</p>

<!-- LICENSE -->
## License
This project is licensed under the [MIT License](https://github.com/TheLuoFengLab/MolLangBench/blob/main/LICENSE). You are free to use, modify, and distribute this codebase under the terms of the MIT license.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
