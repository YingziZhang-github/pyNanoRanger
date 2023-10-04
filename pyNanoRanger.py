#!/usr/bin/env python3
# Yang Liu, Chongwei Bi, Yingzi Zhang, 20231004
# Email: {yang.liu.3, chongwei.bi, yingzi.zhang}@kaust.edu.sa

import os
import re
import sys
import time
import subprocess
import glob
import shutil
import argparse
import logging
import multiprocessing
from _version import version
from datetime import datetime
from timeit import default_timer as timer
import call_variant
import identity_analysis


def get_argparse():
    parser = argparse.ArgumentParser(description='Real-Time analysis of Nanopore data for Covid-19 sequencing.')
    parser.add_argument('-c', '--cfg', type=str, required=True, help='cfg file')
    parser.add_argument('-x', '--cuda', type=str, required=True, help='cuda')
    parser.add_argument('-p', '--path', type=str, required=True, help='path/to/nanopore_result_folder')
    parser.add_argument('-p1', '--primer1', type=str, required=True, help='primer1')
    parser.add_argument('-p2', '--primer2', type=str, required=True, help='primer2')
    parser.add_argument('-s', '--save_path', type=str,
                        help='path/to/saved_folder Default: pyNanoRanger_result in -p PATH folder')
    parser.add_argument('-g', '--guppy_barcoder', type=validate_file,
                        help='Optional: path/to/guppy_barcoder, when offering this parameter, it will do additional '
                             'demultiplexing using guppy_barcoder --require_barcodes_both_ends --trim_barcodes')
    parser.add_argument('-k', '--barcode_kits', type=str,
                        help='barcode kits used, e.g. "EXP-NBD114 EXP-NBD104" it is required when providing -g/--guppy_barcoder')
    
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
    """Check whether `tool` is on PATH and marked as executable."""
    return shutil.which(tool) is not None


def processing_fastq(args):
    #unzip files 
    save_path = args.save_path
    # cmd = """gzip -d {fastq_path}*.gz
    #          """.format(fastq_path=args.path)
    # a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # if a.stdout != b'':
    #     logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    # if a.stderr != b'':
    #     logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))
    
    #step2: merge small fastq files into one big fastq file.
    cmd = """cat {fastq_path}*.fastq \
            > /nanopore2/20211116_1928_MN32433_FAR04274_a0c6de8d_guppy634_sup/example.fastq
            """.format(fastq_path=args.path)
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    
    #step3: extract amplicon.
    #测序的上一步是pcr.pcr是湿实验，需要用两个primer完成。
    #但pcr产生的sequence并不总是完美的两个primer的产物，而是可能会有非特异性的其他奇奇怪怪的sequence.
    #所以这一步是为了提取我们预期的两个primer的产物。
    #对应description里"pyNanoRanger will filter the expected amplicon(s) sequences into analysis."

    # /nanopore2/20211116_1928_MN32433_FAR04274_a0c6de8d_guppy634_sup/example.fastq \
    # --forward ACAGCCTATGCCCCATTTTGG --reverse CGAAGGAGATGGAGGTCGTC \
    # --strict-mode \
    # > example_seqkit_Primer11-5F1.fastq
    cmd = """seqkit amplicon /nanopore2/20211116_1928_MN32433_FAR04274_a0c6de8d_guppy634_sup/example.fastq \
            --forward {primer1} --reverse {primer2} \
            --strict-mode \
            > example_seqkit_Primer11-5F1.fastq
            """.format(primer1=args.primer1, primer2=args.primer2)
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))
    #step4:
    #对应description里"In the analysis, sequences with only one restriction cutting site, 
    #two restriction cutting sites, 
    #and no less than three restriction cutting sites are divided into three folders."

    #step4.1:extract all reads containing cutting sites
    #######在这个示例里，先提取全部含CTGCAG的reads。其他情况下，digestion site可能是其他序列，而不再是CTGCAG
    cmd = """cat \
    example_seqkit_Primer11-5F1.fastq | \
    seqkit grep --threads 40 --by-seq --pattern CTGCAG \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG.fastq
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    #######提取reads的id。后面提取几个fastq文件之间的reads补集时需要用。
    cmd = """grep '^@' example_seqkit_5F1-Primer11.seqkitCTGCAG.fastq | \
    cut -d' ' -f1 | cut -c2- \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG.read_ids.txt
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    #step4.2:extract all reads containing exactly only one cutting site
    #######提取有且仅含1个CTGCAG的reads
    cmd = """cat \
    example_seqkit_5F1-Primer11.seqkitCTGCAG.fastq | \
    seqkit grep --threads 40 --by-seq --use-regexp --pattern "(CTGCAG.*){2}" -v \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG_1hit.fastq
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    cmd = """grep '^@' example_seqkit_5F1-Primer11.seqkitCTGCAG_1hit.fastq | \
    cut -d' ' -f1 | cut -c2- \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG_1hit.read_ids.txt
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    #step4.3:extract all reads containing one or two cutting sites
    #######提取含1或2个CTGCAG的reads

    cmd = """cat \
    example_seqkit_5F1-Primer11.seqkitCTGCAG.fastq | \
    seqkit grep --threads 40 --by-seq --use-regexp --pattern "(CTGCAG.*){3}" -v \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG.1or2hits.fastq
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    cmd = """grep '^@' example_seqkit_5F1-Primer11.seqkitCTGCAG.1or2hits.fastq | \
    cut -d' ' -f1 | cut -c2- \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG.1or2hits.read_ids.txt
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    ##step4.4:通过补集的方法提取有且含2个CTGCAG的reads
    cmd = """grep -vFf example_seqkit_5F1-Primer11.seqkitCTGCAG_1hit.read_ids.txt \
    example_seqkit_5F1-Primer11.seqkitCTGCAG.1or2hits.read_ids.txt \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG.2hits.read_ids.txt
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    cmd = """seqtk subseq \
    example_seqkit_5F1-Primer11.seqkitCTGCAG.1or2hits.fastq \
    example_seqkit_5F1-Primer11.seqkitCTGCAG.2hits.read_ids.txt \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG.2hits.fastq
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    #step4.5:通过补集的方法提取含3个及以上CTGCAG的reads。
    #这部分reads比较复杂，提取后只放在一个单独的folder里供用户参考，不参与后面的数据自动流程。
    cmd = """grep -vFf example_seqkit_5F1-Primer11.seqkitCTGCAG.1or2hits.read_ids.txt \
    example_seqkit_5F1-Primer11.seqkitCTGCAG.read_ids.txt \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG.3orMorehits.read_ids.txt
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))
    
    cmd = """seqtk subseq \
    example_seqkit_5F1-Primer11.seqkitCTGCAG.fastq \
    example_seqkit_5F1-Primer11.seqkitCTGCAG.3orMorehits.read_ids.txt \
    > example_seqkit_5F1-Primer11.seqkitCTGCAG.3orMorehits.fastq
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))
    
    cmd = """canu -d canu_draft -p assm genomeSize=30000 \
    -nanopore-raw example_seqkit_5F1-Primer11.seqkitCTGCAG_1hit.fastq \
    useGrid=False maxThreads=40 minInputCoverage=4 stopOnLowCoverage=4
    """
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))
    
    # cmd = """medaka_consensus -i example_seqkit_5F1-Primer11.seqkitCTGCAG_1hit.fastq \
    # -d ./canu_draft/assm.contigs.fasta \
    # -o consensus -t 40
    # """
    # a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # if a.stdout != b'':
    #     logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    # if a.stderr != b'':
    #     logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

def main():
    args = get_argparse()

    result_folder = args.save_path
    ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
    log_file = result_folder + '/' + ctime + '_rt_nano.log'
    logging.basicConfig(level=logging.DEBUG,
                        format='%(message)s',
                        filename=log_file,
                        filemode='w')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    processing_fastq(args)

    # fastq_regex = args.save_path + '/analyzed_achieve/accumulated_reads/*.fastq'
    # fastq_file_all = glob.glob(str(fastq_regex))

    # if len(fastq_file_all) >= 1:
    #     logging.info("\n%s    Program start ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    #     logging.info("--> Read fastq file from %s/analyzed_achieve/accumulated_reads/*.fastq" % args.save_path)
    #     variant_call(args)

    # else:
    #     sys.stderr.write("ERROR! No fastq file detected in %s/analyzed_achieve/accumulated_reads/\n" %
    #                         args.save_path)
    #     sys.exit(0)

    

if __name__ == '__main__':
    main()

