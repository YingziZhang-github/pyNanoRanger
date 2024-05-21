
# pyNanoRanger

**pyNanoRanger** is a data analysis tool for **NanoRanger**, a sequencing strategy named nanopore-based rapid acquisition of neighboring genomic regions. This tool performs real-time analysis of Nanopore sequencing data for detecting human genomic structural variants. It monitors a specified directory for new FASTQ files, processes them as they are generated, and categorizes reads based on specified primers and cutting sites.

## Features

- Real-time monitoring and processing of Nanopore sequencing data.
- Flexible configuration with user-defined primers and cutting sites.
- Categorizes reads based on the presence and number of cutting sites.
- Logs detailed processing steps and read counts.

## Installation

### Prerequisites

- Python 3.6+
- Required Python packages:
  - argparse
  - logging
  - subprocess
  - datetime
  - glob
  - tabulate
- Required external tools:
  - seqkit
  - guppy

### Steps

1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/pyNanoRanger.git
   cd pyNanoRanger
   ```

2. Install the required Python packages:
   ```sh
   pip install -r requirements.txt
   ```

3. Ensure that `seqkit` and `guppy` are installed and accessible in your system's PATH.

## Usage

### Command-Line Arguments

- `-c`/`--cfg`: Path to the configuration file (required).
- `-x`/`--cuda`: CUDA configuration (required).
- `-p`/`--path`: Path to the folder containing Nanopore sequencing results (required).
- `-p1`/`--primer1`: Sequence of primer 1 (required).
- `-p2`/`--primer2`: Sequence of primer 2 (required).
- `-cs`/`--cutting_site`: Sequence of the cutting site (required).
- `-s`/`--save_path`: Path to the folder where results will be saved (optional, defaults to a subfolder within the specified path).
- `-g`/`--guppy_barcoder`: Path to the guppy barcoder tool for additional demultiplexing (optional).
- `-k`/`--barcode_kits`: Barcode kits used, required if `-g`/`--guppy_barcoder` is provided.

### Example Command

```sh
python pyNanoRanger.py -c config.cfg -x cuda -p /path/to/nanopore_result_folder -p1 ACAGCCTATGCCCCATTTTGG -p2 CGAAGGAGATGGAGGTCGTC -cs CTGCAG -s /path/to/save_folder
```

### Real-Time Output

The script will continuously monitor the specified directory for new FASTQ files and process them in real-time. Example output:

```plaintext
2024-05-21 12:00:00    Program started...
2024-05-21 12:00:00    Monitoring /path/to/nanopore_result_folder for new FASTQ files...

2024-05-21 12:00:10    Detected new FASTQ file: /path/to/nanopore_result_folder/sample_01.fastq

+----------------------------------------------+---------------+---------------------------+-----------------------------+-------------------------------------+-------------------------------------------+
| File                                         | Total Reads   | Reads with Cutting Sites  | Reads with One Cutting Site | Reads with One or Two Cutting Sites | Reads with Three or More Cutting Sites    |
+----------------------------------------------+---------------+---------------------------+-----------------------------+-------------------------------------+-------------------------------------------+
| /path/to/nanopore_result_folder/sample_01.fastq | 15000         | 9000                      | 3000                        | 6000                                | 1000                                       |
+----------------------------------------------+---------------+---------------------------+-----------------------------+-------------------------------------+-------------------------------------------+

Monitoring for new files... (Press Ctrl+C to stop)

2024-05-21 12:01:00    Detected new FASTQ file: /path/to/nanopore_result_folder/sample_02.fastq

+----------------------------------------------+---------------+---------------------------+-----------------------------+-------------------------------------+-------------------------------------------+
| File                                         | Total Reads   | Reads with Cutting Sites  | Reads with One Cutting Site | Reads with One or Two Cutting Sites | Reads with Three or More Cutting Sites    |
+----------------------------------------------+---------------+---------------------------+-----------------------------+-------------------------------------+-------------------------------------------+
| /path/to/nanopore_result_folder/sample_01.fastq | 15000         | 9000                      | 3000                        | 6000                                | 1000                                       |
| /path/to/nanopore_result_folder/sample_02.fastq | 18000         | 12000                     | 4000                        | 8000                                | 2000                                       |
+----------------------------------------------+---------------+---------------------------+-----------------------------+-------------------------------------+-------------------------------------------+

Monitoring for new files... (Press Ctrl+C to stop)

```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or issues, please contact [Yingzi Zhang](mailto:yingzi.zhang@kaust.edu.sa) or [Yang Liu](mailto:yang.liu.3@kaust.edu.sa)
