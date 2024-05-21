# pyNanoRanger
Data analysis tool for NanoRanger.

Usage

command:

python pyNanoRanger.py -c config.cfg -x cuda -p /path/to/nanopore_result_folder -p1 $primer1 -p2 $primer2 -cs $cutting_site -s /path/to/save_folder

options:
-c/--cfg: Path to the configuration file (required).

-x/--cuda: CUDA configuration (required).

-p/--path: Path to the folder containing Nanopore sequencing results (required).

-p1/--primer1: Sequence of primer 1 (required).

-p2/--primer2: Sequence of primer 2 (required).

-cs/--cutting_site: Sequence of the cutting site (required).

-s/--save_path: Path to the folder where results will be saved (optional, defaults to a subfolder within the specified path).

-g/--guppy_barcoder: Path to the guppy barcoder tool for additional demultiplexing (optional).

-k/--barcode_kits: Barcode kits used, required if -g/--guppy_barcoder is provided.
