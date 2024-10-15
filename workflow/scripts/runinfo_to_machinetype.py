#!/usr/bin/env python
from __future__ import print_function  #  support python3-style print 
import sys
import os
import re
import argparse


def get_machine_type(args):
    
    machine_type = "unknown"

    path = args["input_file"][0]
    with open(path,"r") as info:
        for record in info:
            match=re.search("\<Instrument\>([^<]+)\<\/Instrument\>", record)
            if match is None:
                continue
            
            instrument=match.groups()[0]
            
            if instrument[0] == "M":
                machine_type = "miseq"
            elif instrument[0] == "F":
                machine_type = "iseq"
            else:
                machine_type = "novaseq"

            break

    return machine_type

        

def get_options():
    description = """infer machine type from runinfo"""
    long_description = """

examples :

python runinfo_to_machinetype.py /dataset/2023_illumina_sequencing_a/scratch/20230504_FS10001778_51_BSB09412-3119/RunInfo.xml  # iseq
python runinfo_to_machinetype.py /dataset/miseq/active/140130_M02412_0001_000000000-A6PC8/RunInfo.xml                          # miseq
python runinfo_to_machinetype.py /dataset/2023_illumina_sequencing_a/scratch/230619_A01439_0186_AHFHCCDRX3/RunInfo.xml         # novaseq

    """
    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_file', type=str, nargs=1,metavar="path to RunInfo.xml", help='path to RunInfo.xml file')
    parser.add_argument('-t', '--task' , dest='task', required=False, default="get_machine_type" , type=str,
                        choices=["get_machine_type"], help="what you want to do")

    args = vars(parser.parse_args())

    # check files exist
    path = args["input_file"][0]    
    if not os.path.isfile(path) and not os.path.islink(path):
        raise Exception("%s not found, or is not a file or link,  giving up"%path)


    return args



def main():

    opts = get_options()

    if opts["task"] == "get_machine_type":
        print(get_machine_type(opts))
    else:
        raise Exception("oops task %(task)s is not supported yet !"%opts)

if __name__ == "__main__":
   main()

