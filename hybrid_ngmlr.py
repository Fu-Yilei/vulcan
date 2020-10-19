import os
import sys
import subprocess
import argparse
import logging
# import re
import pysam
from tqdm import tqdm
import numpy as np


'''
hybrid_NGMLR: a pipeline connects minimap2 and NGMLR
@author: Yilei Fu
@Email: yf20@rice.edu
'''


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

WORK_DIR = "./work"
THREADS = 1
# LOG = ""


def parseArgs(argv):
    """Function for parsing Arguments

    Args:
        argv: The arguments sent into this program
    Returns:
        arguments
    """

    parser = argparse.ArgumentParser(
        description="hybrid_NGMLR: a pipeline connects minimap2 and NGMLR")

    parser.add_argument("-w", "--WORK_DIR",  type=str,
                        help="Directory of work directory",
                        default="./work")
    parser.add_argument("-t", "--threads",  type=int,
                        help="threads",
                        default=1)
    parser.add_argument("-i", "--input", type=str,
                        help="input read path")
    parser.add_argument("-o", "--output", type=str,
                        help="output bamfile path")
    parser.add_argument("-p", "--percentile", type=int,
                        help="percentile of cut-off, default: 50", default=50)
    parser.add_argument("-r", "--reference", type=str,
                        help="raference path")
    parser.add_argument("-R", "--raw_edit_distance",
                        help="Use raw edit distance to do the cut-off", action="store_true")
    parser.add_argument(
        "-P", "--pacbio", help="Input reads is pacbio reads", action="store_true")
    args = parser.parse_args(argv)

    return args


def run_minimap2(minimap2_input_ref, minimap2_input_reads, minimap2_output,  if_pacbio):
    if if_pacbio:
        minimap2_cmd = f"minimap2 -x map-pb -a -t {THREADS} --MD -o {minimap2_output}  {minimap2_input_ref} {minimap2_input_reads}"
    else:
        minimap2_cmd = f"minimap2 -x map-ont -a -t {THREADS} --MD -o {minimap2_output}  {minimap2_input_ref} {minimap2_input_reads}"

    logger.info(f"Executing: {minimap2_cmd}")

    os.system(minimap2_cmd)
    # subprocess.check_output(minimap2_cmd, shell=True)
    return minimap2_output


def run_ngmlr(ngmlr_input_ref, ngmlr_input_reads, ngmlr_output,  if_pacbio):
    if if_pacbio:
        ngmlr_cmd = f"ngmlr -x pacbio -t {THREADS} -r {ngmlr_input_ref} -q {ngmlr_input_reads} -o {ngmlr_output}"

    else:
        ngmlr_cmd = f"ngmlr -x ont -t {THREADS} -r {ngmlr_input_ref} -q {ngmlr_input_reads} -o {ngmlr_output}"
    logger.info(f"Executing: {ngmlr_cmd}")
    os.system(ngmlr_cmd)
    return ngmlr_output


def keep_primary_mapping(input_sam, output_sam):
    samtools_cmd = f"samtools view -h -@ {THREADS - 1} -F 2308 {input_sam} -o {output_sam}"
    logger.info(f"Executing: {samtools_cmd}")
    os.system(samtools_cmd)
    return output_sam


def generate_distance_file(input_sam, output_txt, raw_edit_distance):
    if raw_edit_distance:
        os.system(f"grep -o -E \"NM:i:[0-9]+\" {input_sam} > {output_txt}")
    else:
        samfile = pysam.AlignmentFile(input_sam, "r")
        with open(output_txt, "w") as out_f:
            for read in tqdm(samfile.fetch()):
                tags = dict(read.tags)
                if "NM" in tags:
                    NM_distance = int(tags["NM"])
                    normalized_edit_distance = float(
                        NM_distance)/read.query_alignment_length
                    out_f.writelines(f"NM:i:{normalized_edit_distance}\n")
    return output_txt


def get_distance_percentile_from_file(distance_file, percentile):
    with open(distance_file, "r") as distance_f:
        distance = distance_f.readlines()
    distance_list = []
    for dist in distance:
        distance_list.append(float(dist[:-1].split(":")[2]))
    # dist_mean = np.mean(distance_list)
    dist_percentile_num = np.percentile(distance_list, percentile)
    return dist_percentile_num


def generate_sorted_bam(input_sam, output_bam):
    transform_cmd = f"samtools view -bS -@ {THREADS} {input_sam} | samtools sort -@ {THREADS} -o {output_bam}"
    logger.info(transform_cmd)
    os.system(transform_cmd)
    return output_bam


def seperate_sam_files(input_sam, under_value_bam, above_value_bam, cut_off_value,  raw_edit_distance):

    transformed_bam = f"{input_sam[:-4]}.bam"
    generate_sorted_bam(input_sam, transformed_bam)
    if raw_edit_distance:
        logger.info(
            f"seperate sam file based on raw edit distance {cut_off_value}")
        bamtools_above_cmd = f"bamtools filter -in {transformed_bam} -out {above_value_bam} -tag NM:{int(cut_off_value)}"
        logger.info(bamtools_above_cmd)
        os.system(bamtools_above_cmd)
        under_filter = "{\"tag\": \"NM:<"+str(cut_off_value)+"\"}"
        under_filter_json = os.path.join(WORK_DIR, "under_filter.json")
        with open(under_filter_json, "w") as ufj:
            ufj.writelines(str(under_filter))
        bamtools_under_cmd = f"bamtools filter -in {transformed_bam} -script {under_filter_json} -out {under_value_bam}"
        logger.info(bamtools_under_cmd)
        os.system(bamtools_under_cmd)
    else:
        logger.info(
            f"Splitting sam file with normalized edit distance {cut_off_value}")
        # above_value_sam = f"{input_sam[:-4]}_above.sam"
        # under_value_sam = f"{input_sam[:-4]}_under.sam"
        samfile = pysam.AlignmentFile(input_sam, "r")
        above_f = pysam.AlignmentFile(above_value_bam, "wb", template=samfile)
        under_f = pysam.AlignmentFile(under_value_bam, "wb", template=samfile)

        # with open(input_sam, "r") as input_f:

        # with open(above_value_sam, "w") as above_f:
        #     with open(under_value_sam, "w") as under_f:
        # for line in input_f:
        #     if line[0] == "@":
        #         above_f.writelines(line)
        #         under_f.writelines(line)
        #     else:
        #         break
        for read in tqdm(samfile.fetch()):
            tags = dict(read.tags)
            if "NM" in tags:
                NM_distance = int(tags["NM"])
                normalized_edit_distance = float(
                    NM_distance)/read.query_alignment_length
                if normalized_edit_distance >= cut_off_value:
                    above_f.write(read)
                else:
                    under_f.write(read)
        # generate_sorted_bam(above_value_sam, above_value_bam)
        # generate_sorted_bam(under_value_sam, under_value_bam)
    # logger.info("sam file split finished")


def bam_to_reads(input_bam, output_reads):
    bam_to_reads_cmd = f"samtools fasta -@ {THREADS} {input_bam} > {output_reads} "
    logger.info(bam_to_reads_cmd)
    os.system(bam_to_reads_cmd)


def merge_bam_files(under_value_bam, above_value_bam, final_output):
    # final_output = os.path.join(WORK_DIR, "hybrid_NGMLR.bam")
    merge_bam_cmd = f"samtools merge {final_output} {under_value_bam} {above_value_bam} -@ {THREADS}"
    logger.info(f"Executing: {merge_bam_cmd}")
    os.system(merge_bam_cmd)

def generate_config_for_notebook():
    #TODO
    return 

def main(argv):
    global THREADS, WORK_DIR
    args = parseArgs(argv)
    WORK_DIR = args.WORK_DIR
    try:
        os.mkdir(WORK_DIR)
    except Exception:
        logging.exception("work path already exists!")
    THREADS = args.threads
    read_file = args.input
    reference_file = args.reference
    percentile = args.percentile
    raw_edit_distance = args.raw_edit_distance

    pacbio = args.pacbio
    final_output = args.output
    print(f"RAW: {raw_edit_distance}, pacbio: {pacbio}")
    # LOG = os.path.join(WORK_DIR, "log.log")
    minimap2_full_sam = os.path.join(WORK_DIR, "minimap2_full.sam")
    minimap2_full_sam_primary = os.path.join(
        WORK_DIR, "minimap2_full_primary.sam")
    minimap2_full_distance = os.path.join(
        WORK_DIR, "minimap2_full_distance.txt")
    # minimap2_under_percentile_sam = os.path.join(
    #     WORK_DIR, f"minimap2_under{percentile}.sam")
    minimap2_above_percentile_bam = os.path.join(
        WORK_DIR, f"minimap2_above{percentile}.bam")
    minimap2_under_percentile_bam = os.path.join(
        WORK_DIR, f"minimap2_under{percentile}.bam")
    above_percentile_reads = os.path.join(
        WORK_DIR, f"above_{percentile}.fasta")
    ngmlr_output = os.path.join(
        WORK_DIR, f"ngmlr_above{percentile}.bam")

    logger.info("run minimap2 on entire reads")
    run_minimap2(reference_file, read_file, minimap2_full_sam, pacbio)
    logger.info("...finished")
    logger.info("keep the primary mapping for edit distance calculation")
    keep_primary_mapping(minimap2_full_sam, minimap2_full_sam_primary)
    logger.info("...finished")
    logger.info("generate edit distance file")
    generate_distance_file(minimap2_full_sam_primary,
                           minimap2_full_distance, raw_edit_distance)
    logger.info("...finished")
    logger.info("gettng the number of cut-off")
    dist_percentile_num = get_distance_percentile_from_file(
        minimap2_full_distance, percentile)
    logger.info("...finished")
    logger.info("seperate sam files")
    seperate_sam_files(minimap2_full_sam, minimap2_under_percentile_bam,
                       minimap2_above_percentile_bam, dist_percentile_num, raw_edit_distance)
    logger.info("...finished")
    logger.info("extracting reads have edit distance above cut-off value")
    bam_to_reads(minimap2_above_percentile_bam,
                 above_percentile_reads)
    logger.info("...finished")
    logger.info("run NGMLR on extracted reads")
    run_ngmlr(reference_file, above_percentile_reads,
              ngmlr_output, pacbio)
    logger.info("...finished")
    logger.info("merge NGMLR's result and minimap2's under cut-off results")
    merge_bam_files(minimap2_under_percentile_bam, ngmlr_output, final_output)
    logger.info("...finished")
    logger.info("All Finished")


if __name__ == "__main__":
    main(sys.argv[1:])
