# COVID Genome Variation

#### Tenderini Ruben - 879290
#### Falbo Andrea - 887525

## Table of Contents

#### 1. [Objective](#objective)
#### 2. [Contents](#contents)
#### 3. [Tools](#tools)
#### 4. [Setup](#setup)
#### 5. [Codebase](#codebase)
#### 6. [Usage](#usage)
#### 6.  [Collaborate](#collaborate)

## Objective
This project involves a Python script that utilizes MAFFT to align SARS-CoV-2 genomes concerning the reference genome, 
identifying and reporting pointwise variations compared to the reference, along with comparative statistics on the analyzed genomes.

## Contents

The goal is to develop a Python script capable of analyzing the input file `covid-sequences.fasta`, comprising SARS-CoV-2 genomes sequenced in November 2021 and retrieved from the NCBI website. The initial genome, identified as NC_045512.2, acts as the reference genome sequenced in autumn 2019.

We employ MAFFT for alignment, leveraging the resulting multiple sequence alignment matrix to identify all pointwise variations among genomes in comparison to the reference. The script is designed to generate a comprehensive report detailing these variations, including their types (substitution, insertion in the reference, deletion in the reference), the bases involved (or the inserted/deleted base), and the count of genomes exhibiting the variation relative to the reference. Notably, substitutions involving the 'N' base are disregarded, whereas insertions/deletions involving 'N' are accounted for. Additionally, the report highlights:
- The genome exhibiting the most and least variations in comparison to the reference.
- Positions in the reference where variations occur across different genomes.
- Positions in the reference where variations occur identically across different genomes.


## Tools

1. **MAFFT** stands for Multiple Alignment using Fast Fourier Transform. It's a bioinformatics tool used for the 
alignment of multiple DNA, RNA, or protein sequences. MAFFT is commonly used to align sequences in comparative genomics 
studies to understand evolutionary relationships or identify variations between sequences.

2. **Python** is a high-level, interpreted programming language known for its simplicity and readability. It is widely used in bioinformatics and computational biology due to its versatility and extensive libraries like Biopython, NumPy, and SciPy. Python allows bioinformaticians to develop powerful and efficient tools for analyzing biological data and solving complex problems.

3. **PyCharm** is an integrated development environment (IDE) specifically designed for Python development. It provides features such as code debugging, syntax highlighting, intelligent code completion, and project management tools. PyCharm helps bioinformaticians and developers write, test, and debug Python scripts more efficiently.

4. **GitHub** is a web-based platform and version control system used for hosting and sharing code repositories. It allows developers to collaborate on projects, track changes, manage issues, and merge code changes through pull requests. GitHub is widely used in the bioinformatics community for sharing bioinformatics tools, pipelines, and data analysis workflows. It fosters collaboration and facilitates the dissemination of scientific research and software development in bioinformatics.

## Setup

1. **Install Python:** Download and install Python from the official website: [Python.org](https://www.python.org/).

If you're using PyCharm, after following these simple steps, everything should be ready to go!

The same should apply for VSCode!

## Codebase
- *docs* directory containing:
  - `Progetto3.pdf`: the project assignment file.
- *fasta* directory containing input fasta files:
  - `covid-sequences.fasta`: the file assigned in the project description.
  - `aligned-covid-sequences.fasta`: fasta files with aligned strings via MAFFT.
- *report* directory containing output files generated:
  - `comparisons.txt`: contains comparisons between sequences.
  - `final_report.txt`: contains the final requested results.
- *scripts* containing necessary Python script for the program:
  - `analyzer.py`: accepts a FASTA file with genome sequences and:
    - validates sequences for length and base composition
    - identifies variations compared to a reference sequence 
    - generates a final report with variation details
- `main.py`: demonstrates the usage of the GenomeVariationAnalyzer class defined in `analyzer.py`. 

## Usage

If you wish to use this project with your custom files, follow these steps:

1. **Prepare Custom Files:**
   - Ensure you have a FASTA file containing the genome sequences you want to analyze.
   
2. **Replace Default FASTA File:**
   - Replace the default FASTA file path (`fasta/aligned-covid-sequences.fasta`) in the `main.py` file with the path to your custom FASTA file.

3. **Run the Program:**
   - Once you've replaced the default FASTA file with your custom file, run the `main.py` file.
   - This will initiate the analysis process using your custom FASTA file and generate the corresponding output reports.

4. **Explore the Results:**
   - After the program execution completes, you can examine the generated output reports in the `report` directory.
   - The output files will include `comparisons.txt`, which displays differences between the analyzed sequences, and `final_report.txt`, containing the final analysis results.

By following these steps, you'll be able to use the project with your custom files and analyze the genome sequences of your interest. If you encounter any issues or have questions, feel free to reach out to the development team for assistance and support.


## Collaborate
If you encounter any errors or have suggestions for improvements, we welcome collaboration and feedback from the community. You can contribute by:
- **Reporting Issues:** If you come across any bugs or issues, please submit them through GitHub's issue tracker for this project repository.
- **Pull Requests:** Feel free to submit pull requests with fixes, enhancements, or new features. We appreciate any contributions that improve the project.
Collaboration is essential for the continued development and improvement of this project. Let's work together to make it even better!