
import os
import sys
import time
import subprocess
import argparse
import logging
import glob
from datetime import datetime
from _version import version
from tabulate import tabulate

import call_variant
import identity_analysis


def get_argparse():
    parser = argparse.ArgumentParser(description='Real-Time analysis of Nanopore data for detecting human genomic structural variants.')
    parser.add_argument('-c', '--cfg', type=str, required=True, help='cfg file')
    parser.add_argument('-x', '--cuda', type=str, required=True, help='cuda')
    parser.add_argument('-p', '--path', type=str, required=True, help='path/to/nanopore_result_folder')
    parser.add_argument('-p1', '--primer1', type=str, required=True, help='primer1')
    parser.add_argument('-p2', '--primer2', type=str, required=True, help='primer2')
    parser.add_argument('-cs', '--cutting_site', type=str, required=True, help='cutting site sequence')
    parser.add_argument('-s', '--save_path', type=str, help='path/to/saved_folder Default: pyNanoRanger_result in -p PATH folder')
    parser.add_argument('-g', '--guppy_barcoder', type=validate_file, help='Optional: path/to/guppy_barcoder, when offering this parameter, it will do additional demultiplexing using guppy_barcoder --require_barcodes_both_ends --trim_barcodes')
    parser.add_argument('-k', '--barcode_kits', type=str, help='barcode kits used, required if -g/--guppy_barcoder is provided')
    
    args = parser.parse_args()
    if args.guppy_barcoder is not None:
        if args.barcode_kits is None:
            sys.stderr.write("ERROR!  Please provide -k/--barcode_kits\n")
            sys.exit(1)
        else:
            kit_list = barcode_kit_list()
            for kit in str(args.barcode_kits).split(" "):
                if kit not in kit_list:
                    print_kit_list = "\n".join(kit_list)
                    sys.stderr.write("ERROR!  -k/--barcode_kits is not in kit list:\n\n%s\n\n" % print_kit_list)
                    sys.exit(1)
    return args


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def barcode_kit_list():
    choices = ["EXP-NBD103", "EXP-NBD104", "EXP-NBD114", "EXP-NBD196", "EXP-PBC001", "EXP-PBC096",
               "OND-SQK-LP0096M", "OND-SQK-LP0096S", "OND-SQK-LP1152S", "OND-SQK-LP9216",
               "SQK-16S024", "SQK-LWB001", "SQK-PBK004", "SQK-PCB109", "SQK-RAB201", "SQK-RAB204",
               "SQK-RBK001", "SQK-RBK004", "SQK-RBK096", "SQK-RLB001", "SQK-RPB004",
               "VSK-VMK001", "VSK-VMK002"]
    return choices


def check_tool(tool):
    return shutil.which(tool) is not None


def count_reads(file_path):
    cmd = f"grep -c '^@' {file_path}"
    result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return int(result.stdout.strip())


def process_fastq_file(file_path, save_path, primer1, primer2, cutting_site, summary):
    try:
        merged_fastq = os.path.join(save_path, "merged.fastq")
        cmd = f"cat {file_path}/*.fastq > {merged_fastq}"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        extracted_fastq = os.path.join(save_path, "extracted.fastq")
        cmd = f"seqkit amplicon {merged_fastq} --forward {primer1} --reverse {primer2} --strict-mode > {extracted_fastq}"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        total_reads = count_reads(extracted_fastq)
        
        contains_cutting_sites_fastq = os.path.join(save_path, "contains_cutting_sites.fastq")
        cmd = f"cat {extracted_fastq} | seqkit grep --threads 40 --by-seq --pattern {cutting_site} > {contains_cutting_sites_fastq}"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        reads_with_cutting_sites = count_reads(contains_cutting_sites_fastq)
        
        one_cutting_site_fastq = os.path.join(save_path, "one_cutting_site.fastq")
        cmd = f"cat {contains_cutting_sites_fastq} | seqkit grep --threads 40 --by-seq --use-regexp --pattern '({cutting_site}.*){{2}}' -v > {one_cutting_site_fastq}"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        reads_one_cutting_site = count_reads(one_cutting_site_fastq)

        one_or_two_cutting_sites_fastq = os.path.join(save_path, "one_or_two_cutting_sites.fastq")
        cmd = f"cat {contains_cutting_sites_fastq} | seqkit grep --threads 40 --by-seq --use-regexp --pattern '({cutting_site}.*){{3}}' -v > {one_or_two_cutting_sites_fastq}"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        reads_one_or_two_cutting_sites = count_reads(one_or_two_cutting_sites_fastq)

        three_or_more_cutting_sites_fastq = os.path.join(save_path, "three_or_more_cutting_sites.fastq")
        cmd = f"grep -vFf {one_or_two_cutting_sites_fastq} {contains_cutting_sites_fastq} > {three_or_more_cutting_sites_fastq}"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        reads_three_or_more_cutting_sites = count_reads(three_or_more_cutting_sites_fastq)
        
        summary.append([file_path, total_reads, reads_with_cutting_sites, reads_one_cutting_site, reads_one_or_two_cutting_sites, reads_three_or_more_cutting_sites])
        print_summary(summary)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error processing file {file_path}: {e}")


def print_summary(summary):
    headers = ["File", "Total Reads", "Reads with Cutting Sites", "Reads with One Cutting Site", "Reads with One or Two Cutting Sites", "Reads with Three or More Cutting Sites"]
    print("\n" + tabulate(summary, headers=headers, tablefmt="grid") + "\n")


def continuous_processing(args):
    processed_files = set()
    summary = []
    while True:
        fastq_files = glob.glob(os.path.join(args.path, '*.fastq'))
        new_files = [f for f in fastq_files if f not in processed_files]
        
        if new_files:
            for file_path in new_files:
                process_fastq_file(file_path, args.save_path, args.primer1, args.primer2, args.cutting_site, summary)
                processed_files.add(file_path)
        
        print("Monitoring for new files... (Press Ctrl+C to stop)")
        time.sleep(10)  # Check for new files every 10 seconds


def main():
    args = get_argparse()

    result_folder = args.save_path or os.path.join(args.path, "pyNanoRanger_result")
    os.makedirs(result_folder, exist_ok=True)

    log_file = os.path.join(result_folder, datetime.now().strftime("%Y%m%d_%H.%M.%S") + '_rt_nano.log')
    logging.basicConfig(level=logging.DEBUG, format='%(message)s', filename=log_file, filemode='w')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    continuous_processing(args)


if __name__ == '__main__':
    main()
