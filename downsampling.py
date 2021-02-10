#!/usr/bin/env python3

import logging
import argparse
import csv
import re
import subprocess
import os
import sys
import pathlib
import glob

def parse_args():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True, default=None,
                        help="Input CRAM or BAM file. Default: None.")
    parser.add_argument("--output_dir", required=False, default="./",
                        help="Path to output directory. Default: ./")
    parser.add_argument("--coverage", required=False, default="15",
                        help="Desired coverage for new BAM file. Default: 15") 
    parser.add_argument("--thread", required=False, default="4",
                        help="Number of threads to use (at least 4 for Mosdepth...deal with it). Default: 4")                    
    parser.add_argument("--sample_name", required=False, default=None,
                        help="Sample Name to use, otherwise will default to CRAM files sample name. Useful when dealing with benchmarking data to avoid overwriting eg. SRR11321732_hg19 vs SRR11321732_GRCh38.")
    parser.add_argument("--skip_dupRm", required=False, default=False, action="store_true",
                        help="Do not do a dup remove (directly downsample original BAM). Default: False.")
    parser.add_argument("--tmp", required=False, default=None,
                        help="Temporary directory, deleted after processing. Default: ./sampleID/tmp ")
    parser.add_argument("--loglevel", required=False, default="DEBUG",
                        help="Set logging level to DEBUG (default), INFO or WARNING.")
    args = parser.parse_args()

    # checks
    if not os.path.exists(args.input_file):
        logging.error("Input file not found: %s" % args.input_file)
        sys.exit(1)

    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    try:
        subprocess.call(["which", "mosdepth"], stdout=subprocess.DEVNULL)
    except:
        logging.error("mosdepth not found")
        sys.exit(1)

    try:
        subprocess.call(["which", "sambamba"], stdout=subprocess.DEVNULL)
    except:
        logging.error("sambamba not found")
        sys.exit(1)

    return args
    
def set_logging(loglevel):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % loglevel)
    logging.basicConfig(format="%(asctime)s %(levelname)-8s %(message)s", level=numeric_level)

def check_type(input_file):
    if input_file.endswith(('.cram','.CRAM')):
        return "cram"
    elif input_file.endswith(('.bam','.BAM')):
        return "bam"
    else:
        # Not sure what are you processing, must DIE
        logging.error("input_file is not bam or cram. Dying now!")
        sys.exit(1)

def get_ID(input_file , sample_name):
    if sample_name is not None:
        return sample_name
    else:
        # assume sample ID is the first segment (sample1) in bam file named sample1.bqsr.sorted.whatever.bam
        bam = input_file.split('/')[-1]
        sampleID = bam.split('.')[0]
        return sampleID

def temp_output_loc( out_dir , tmp , sampleID , cov):
    if tmp is not None:
        tempLoc = "%s/%s/%sx" % (tmp , sampleID , cov)
        logging.info("Creating temp folder %s" % tempLoc)
        pathlib.Path(tempLoc).mkdir(parents=True, exist_ok=True)
        return tempLoc
    else:
        tempLoc = "%s/%s/tmp" % (out_dir , sampleID)
        logging.info("Creating temp folder %s" % tempLoc)
        pathlib.Path(args.out_dir).mkdir(parents=True, exist_ok=True)
        return tempLoc

def remove_dup(input_file , tempLoc , file_type , sampleID , skip_dupRm , thread):
    filterString = "not (unmapped or mate_is_unmapped) and paired and not duplicate"
    cram = "-C"
    logging.info("Sanitising input %s " % input_file)
    if skip_dupRm :
        logging.info("Not removing duplicates reads")
        filterString = "not (unmapped or mate_is_unmapped) and paired"       
    if re.match("bam", file_type):
        cram = ""
        
    # run 
    # sambamba view -F 'not (unmapped or mate_is_unmapped) and paired and not duplicate' -q -C -f bam -t ${THREAD} -o ${OUT}/${SAMPLE}.dupRm.bam ${INPUT_CRAM}
    cmd = "sambamba view -F '%s' -q %s -f bam -t %d -o %s/%s.dupRm.bam %s" % ( filterString , cram , int(thread) , tempLoc , sampleID , input_file )
    if os.system(cmd) != 0:
        logging.error("Sambamba sanitising bam/cram failed. Dying now.")
        sys.exit(1)
    else:
        bamLoc = "%s/%s.dupRm.bam" % ( tempLoc , sampleID )
        logging.info("Sanitised bam found at %s " % bamLoc)
        return bamLoc

def get_coverage( Bam , sampleID , tempFolder , subfolder ):
    # run
    # mosdepth --no-per-base --threads 4 ${OUT_L} ${IN}
    sub_tmp = "%s/%s" % ( tempFolder , subfolder )
    pathlib.Path(sub_tmp).mkdir(parents=True, exist_ok=True)
    cmd = "mosdepth --no-per-base --threads 4 %s/%s %s " % ( sub_tmp , sampleID , Bam) 
    if os.system(cmd) != 0:
        logging.error("Get coverage %s failed. Dying now." % subfolder )
        sys.exit(1)
    
    summary_file = "%s/%s.mosdepth.summary.txt" % ( sub_tmp , sampleID)
    mos_cov = open(summary_file, "r")
    length = 0
    bases = 0 
    for line in mos_cov :
        # only counts for chr1-22 or 1-22
        if re.search(r"^[0-9]{1,2}", line) or re.search(r"^chr[0-9]{1,2}", line):
            array_line = re.split(r'\t+', line.strip())
            length += int(array_line[1])
            bases += int(array_line[2])
    
    cov = bases/length
    logging.info("%s coverage is : %f" % ( Bam , cov ) )
    return cov

def get_fraction( neededCov , bamCov):
    if neededCov >= bamCov :
        logging.error("Needed coverage ( %f ) is more than the CRAM's coverage ( %f ): unable to downsample." % ( neededCov , bamCov))
        sys.exit(1)
    else:
        fraction = float( neededCov ) / float( bamCov ) 
        logging.info("Downsampling CRAM with factor: %f " % fraction)
        return fraction

def downsample_bam(dupRemovedBam , output_dir , sampleID , fraction , seed , cov ):
    # run
    # sambamba view -s ${FACTOR} --subsampling-seed=${SEED} -f bam -o ${OUT_DOWNSAMPLING}/${SAMPLE}/${SAMPLE}.dupRm.subsam.bam ${OUT}/${SAMPLE}.dupRm.bam
    out_prefix = "%s/%s/%sx" % ( output_dir , sampleID , cov )
    pathlib.Path(out_prefix).mkdir(parents=True, exist_ok=True)
    cmd = "sambamba view -s %f --subsampling-seed=%d -f bam -o %s/%s.dupRm.subsam.bam %s" % ( fraction , seed , out_prefix , sampleID , dupRemovedBam)
    logging.info("Downsampling %s " % dupRemovedBam )
    if os.system(cmd) != 0:
        logging.error("Sambamba downsampling failed. Dying now.")
        sys.exit(1)
    else:
        bamLoc = "%s/%s.dupRm.subsam.bam" % ( out_prefix , sampleID )
        logging.info("Downsampled bam found at %s " % bamLoc)
        return bamLoc

def cleanUp ( tempFolder ):
    logging.info("Cleaning up, deleting %s" % tempFolder)
    cmd = "rm -rf %s" % tempFolder
    os.system(cmd)
    logging.info("Done!")

if __name__ == "__main__":
    # seed number to make result reproducible
    seed = 10 
    
    args = parse_args()
    set_logging(args.loglevel)

    # check input is CRAM or BAM
    file_type = check_type(args.input_file)
    # get sample ID
    sampleID = get_ID(args.input_file  , args.sample_name )
    # check and create intermediate ouput directory
    tempFolder = temp_output_loc ( args.output_dir , args.tmp , sampleID , args.coverage)
    # Remove dup and keep ONLY paired reads + inflat if CRAM
    dupRemovedBam = remove_dup(args.input_file , tempFolder , file_type , sampleID , args.skip_dupRm , args.thread)
    # Get coverage of input file - Round 1
    cov = get_coverage( dupRemovedBam , sampleID , tempFolder , "R1")
    # Check Coverage of input is larger than desired coverage, if true get "fraction" needed for downsampling
    fraction = get_fraction( args.coverage , cov )
    # Downsample to fraction
    downsampledBam = downsample_bam(dupRemovedBam , args.output_dir , sampleID , fraction , seed , args.coverage )
    # Get coverage for NEW BAM
    newCov = get_coverage( downsampledBam , sampleID , tempFolder , "R2")
    # clean up and Done
    cleanUp( tempFolder )