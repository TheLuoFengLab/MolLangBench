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

## Table of Contents
- [About The Project](#about-the-project)
- [Getting Started](#getting-started)
- [Dataset Structure](#dataset-structure)
- [Usage](#usage)
- [Evaluation](#evaluation)
- [Citation](#citation)
- [Contact](#contact)
- [License](#license)

## About The Project

**MolLangBench** is a comprehensive benchmark designed to evaluate the fundamental capabilities of AI models in language-prompted molecular structure recognition, editing, and generation.

This repository provides:
- Code and examples to load and use the dataset directly from the [Hugging Face Dataset](https://huggingface.co/datasets/ChemFM/MolLangBench)
- Evaluation scripts and templates to test OpenAI models (e.g., o1, o3, o4-mini) using either molecular images or SMILES strings as inputs

It is straightforward to extend this repository to evaluate other language or multimodal models by adapting the provided input formatting and evaluation templates.

## Getting Started

You can easily set up the required environment using Conda by following these steps:

* Clone the repository

    ```bash
    git clone https://github.com/TheLuoFengLab/MolLangBench.git
    cd MolLangBench
    ```

* Create and activate Conda environment
    ```bash
    conda env create -f environment.yml 
    conda activate MolLangBench
    ```
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

### Quick Start

Once the environment is set up, you can load the MolLangBench dataset directly from the Hugging Face Hub:

```python
from datasets import load_dataset

# Load the recognition task (train + test splits available)
rec_train = load_dataset("ChemFM/MolLangBench", name="recognition", split="train")
rec_test = load_dataset("ChemFM/MolLangBench", name="recognition", split="test")

# Load the editing task (only test split)
edit = load_dataset("ChemFM/MolLangBench", name="edit", split="test")

# Load the generation task (only test split)
gen = load_dataset("ChemFM/MolLangBench", name="generation", split="test")
```

<!-- LICENSE -->
## License
This project is licensed under the [MIT License](https://github.com/TheLuoFengLab/MolLangBench/blob/main/LICENSE). You are free to use, modify, and distribute this codebase under the terms of the MIT license.

<p align="right">(<a href="#readme-top">back to top</a>)</p>