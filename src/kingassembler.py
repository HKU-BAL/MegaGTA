#!/usr/bin/env python

from __future__ import print_function

import sys
import getopt
import subprocess
import errno
import os, glob
import shutil
import locale
import signal
import multiprocessing
import logging
import time
from datetime import datetime

# if len(sys.argv) < 4:
#     sys.exit('kingassembler <k-size> <reads-file> <output-path>')
# print 'Argument List:', str(sys.argv)
# if not os.path.exists(sys.argv[3]):
#     os.makedirs(sys.argv[3])
# call(["/nas5/ykhuang/megahit/megahit --graph-only --k-list ", sys.argv[1], " -r ", sys.argv[2], " -o ", sys.argv[3], " --min-count 1"])
# # call(["../megahit/megahit"])

megahit_version_str = ""
usage_message = '''
Copyright (c) The University of Hong Kong

Usage:
  megahit [options] {-1 <pe1> -2 <pe2> | --12 <pe12> | -r <se>} [-o <out_dir>]

  Input options that can be specified for multiple times (supporting plain text and gz/bz2 extensions)
    -1                       <pe1>          comma-separated list of fasta/q paired-end #1 files, paired with files in <pe2>
    -2                       <pe2>          comma-separated list of fasta/q paired-end #2 files, paired with files in <pe1>
    --12                     <pe12>         comma-separated list of interleaved fasta/q paired-end files
    -r/--read                <se>           comma-separated list of fasta/q single-end files

  Input options that can be specified for at most ONE time (not recommended):
    --input-cmd              <cmd>          command that outputs fasta/q reads to stdout; taken by MEGAHIT as SE reads 

Optional Arguments:
  Basic assembly options:
    --min-count              <int>          minimum multiplicity for filtering (k_min+1)-mers, default 2
    --k-min                  <int>          minimum kmer size (<= 127), must be odd number, default 21
    --k-max                  <int>          maximum kmer size (<= 127), must be odd number, default 99
    --k-step                 <int>          increment of kmer size of each iteration (<= 28), must be even number, default 10
    --k-list                 <int,int,..>   comma-separated list of kmer size (all must be odd, in the range 15-127, increment <= 28);
                                            override `--k-min', `--k-max' and `--k-step'

  Advanced assembly options:
    --no-mercy                              do not add mercy kmers
    --no-bubble                             do not merge bubbles
    --merge-level            <l,s>          merge complex bubbles of length <= l*kmer_size and similarity >= s, default 20,0.98
    --prune-level            <int>          strength of local low depth pruning (0-2), default 2
    --low-local-ratio        <float>        ratio threshold to define low local coverage contigs, default 0.2
    --max-tip-len            <int>          remove tips less than this value; default 2*k for iteration of kmer_size=k
    --no-local                              disable local assembly
    --kmin-1pass                            use 1pass mode to build SdBG of k_min

  Presets parameters:
    --presets                <str>          override a group of parameters; possible values:
                                            meta            '--min-count 2 --k-list 21,41,61,81,99'             (generic metagenomes, default)
                                            meta-sensitive  '--min-count 2 --k-list 21,31,41,51,61,71,81,91,99' (more sensitive but slower)
                                            meta-large      '--min-count 2 --k-list 27,37,47,57,67,77,87'       (large & complex metagenomes, like soil)
                                            bulk            '--min-count 3 --k-list 31,51,71,91,99 --no-mercy'  (experimental, standard bulk sequencing with >= 30x depth)
                                            single-cell     '--min-count 3 --k-list 21,33,55,77,99,121 --merge_level 20,0.96' (experimental, single cell data)

  Hardware options:
    -m/--memory              <float>        max memory in byte to be used in SdBG construction; default 0.9
                                            (if set between 0-1, fraction of the machine's total memory)
    --mem-flag               <int>          SdBG builder memory mode, default 1
                                            0: minimum; 1: moderate; others: use all memory specified by '-m/--memory'.
    --use-gpu                               use GPU
    --gpu-mem                <float>        GPU memory in byte to be used. Default: auto detect to use up all free GPU memory. 
    -t/--num-cpu-threads     <int>          number of CPU threads, at least 2. Default: auto detect to use all CPU threads.

  Output options:
    -o/--out-dir             <string>       output directory, default ./megahit_out
    --out-prefix             <string>       output prefix (the contig file will be OUT_DIR/OUT_PREFIX.contigs.fa)
    --min-contig-len         <int>          minimum length of contigs to output, default 200
    --keep-tmp-files                        keep all temporary files

Other Arguments:
    --continue                              continue a MEGAHIT run from its last available check point.
                                            please set the output directory correctly when using this option.
    -h/--help                               print the usage message
    -v/--version                            print version
    --verbose                               verbose mode
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Options():
    def __init__(self):
        self.host_mem = 0.9
        self.gpu_mem = 0
        self.out_dir = "./megahit_out/"
        self.min_contig_len = 200
        self.k_min = 21
        self.k_max = 99
        self.k_step = 20
        self.k_list = list()
        self.k_current = -1
        self.set_list_by_min_max_step = True
        self.min_count = 2
        self.bin_dir = sys.path[0] + "/"
        self.max_tip_len = -1
        self.no_mercy = False
        self.no_local = False
        self.no_bubble = False
        self.merge_len = 20
        self.merge_similar = 0.98
        self.prune_level = 2
        self.num_cpu_threads = False
        self.low_local_ratio = 0.2
        self.temp_dir = self.out_dir + "tmp/"
        self.contig_dir = self.out_dir + "intermediate_contigs/"
        self.keep_tmp_files = False
        self.use_gpu = False
        self.mem_flag = 1
        self.continue_mode = False;
        self.last_cp = -1;
        self.out_prefix = ""
        self.kmin_1pass = False
        self.verbose = False
        self.pe1 = []
        self.pe2 = []
        self.pe12 = []
        self.aligned_ref = []
        self.se = []
        self.input_cmd = ""
        self.inpipe = dict()
        self.presets = ""
        self.graph_only = False
        self.combined_contigs_file = ""
        self.filtered_nucl_file = "nucl_merged.fasta"
        self.filtered_prot_file = "prot_merged.fasta"
        self.filter_len = 450
        self.clustering_java_heap_memory = 16
        self.clustering = "Clustering.jar"
        self.gene_info = {}

opt = Options()
cp = 0

def log_file_name():
    if opt.out_prefix == "":
        return opt.out_dir + "log"
    else:
        return opt.out_dir + opt.out_prefix + ".log"

def opt_file_name():
    return opt.out_dir + "opts.txt"

def make_out_dir():
    if os.path.exists(opt.out_dir):
        pass
    else:
        os.mkdir(opt.out_dir)

    if os.path.exists(opt.temp_dir):
        pass
    else:
        os.mkdir(opt.temp_dir)

    if os.path.exists(opt.contig_dir):
        pass
    else:
        os.mkdir(opt.contig_dir)

def parse_opt(argv):
    try:
        opts, args = getopt.getopt(argv, "hm:o:r:t:v1:2:l:", 
                                    ["help",
                                     "read=",
                                     "12=",
                                     "input-cmd=",
                                     "memory=",
                                     "out-dir=",
                                     "min-contig-len=",
                                     "use-gpu",
                                     "num-cpu-threads=",
                                     "gpu-mem=",
                                     "kmin-1pass",
                                     "k-min=",
                                     "k-max=",
                                     "k-step=",
                                     "k-list=",
                                     "num-cpu-threads=",
                                     "min-count=",
                                     "no-mercy",
                                     "no-local",
                                     "max-tip-len=",
                                     "no-bubble",
                                     "prune-level=",
                                     "merge-level=",
                                     "low-local-ratio=",
                                     "keep-tmp-files",
                                     "mem-flag=",
                                     "continue",
                                     "version",
                                     "out-prefix=",
                                     "verbose",
                                     "presets=",
                                     "graph-only",
                                     "max-read-len=",
                                     "no-low-local",
                                     "cpu-only",
                                     "forward_hmm=",
                                     "reverse_hmm=",
                                     "gene_list="])
    except getopt.error as msg:
        raise Usage(megahit_version_str + '\n' + str(msg))
    if len(opts) == 0:
        raise Usage(megahit_version_str + '\n' + usage_message)

    global opt
    need_continue = False

    for option, value in opts:
        if option in ("-h", "--help"):
            print(megahit_version_str + '\n' + usage_message)
            exit(0)
        elif option in ("-o", "--out-dir"):
            if opt.continue_mode == 0:
                opt.out_dir = value + "/"
        elif option in ("-m", "--memory"):
            opt.host_mem = float(value)
        elif option == "--gpu-mem":
            opt.gpu_mem = long(float(value))
        elif option == "--min-contig-len":
            opt.min_contig_len = int(value)
        elif option in ("-t", "--num-cpu-threads"):
            opt.num_cpu_threads = int(value)
        elif option == "--kmin-1pass":
            opt.kmin_1pass = True
        elif option == "--k-min":
            opt.k_min = int(value)
        elif option == "--k-max":
            opt.k_max = int(value)
        elif option == "--k-step":
            opt.k_step = int(value)
        elif option == "--k-list":
            opt.k_list = list(map(int, value.split(",")))
            opt.k_list.sort()
            opt.set_list_by_min_max_step = False
        elif option == "--min-count":
            opt.min_count = int(value)
        elif option == "--max-tip-len":
            opt.max_tip_len = int(value)
        elif option == "--merge-level":
            (opt.merge_len, opt.merge_similar) = map(float, value.split(","))
            opt.merge_len = int(opt.merge_len)
        elif option == "--prune-level":
            opt.prune_level = int(value)
        elif option == "--no-bubble":
            opt.no_bubble = True
        elif option == "--no-mercy":
            opt.no_mercy = True
        elif option == "--no-local":
            opt.no_local = True
        elif option == "--low-local-ratio":
            opt.low_local_ratio = float(value)
        elif option == "--keep-tmp-files":
            opt.keep_tmp_files = True
        elif option == "--mem-flag":
            opt.mem_flag = int(value)
        elif option in ("-v", "--version"):
            print(megahit_version_str)
            exit(0)
        elif option == "--verbose":
            opt.verbose = True
        elif option == "--continue":
            if opt.continue_mode == 0: # avoid check again again again...
                need_continue = True
        elif option == "--out-prefix":
            opt.out_prefix = value
        elif option in ("--cpu-only", "-l", "--max-read-len", "--no-low-local"):
            continue # historical options, just ignore
        elif option in ("-r", "--read"):
            opt.se += value.split(",")
        elif option == "--forward_hmm":
            opt.for_hmm = value
        elif option == "--reverse_hmm":
            opt.rev_hmm = value
        elif option == "-1":
            opt.pe1 += value.split(",")
        elif option == "-2":
            opt.pe2 += value.split(",")
        elif option == "--12":
            opt.pe12 += value.split(",")
        elif option == "--input-cmd":
            opt.input_cmd = value
        elif option == "--presets":
            opt.presets = value
        elif option == "--graph-only":
            opt.graph_only = True;
        elif option == "--gene_list":
            opt.gene_list = value

        else:
            raise Usage("Invalid option %s", option)

    opt.temp_dir = opt.out_dir + "tmp/"
    opt.contig_dir = opt.out_dir + "intermediate_contigs/"

    if need_continue:
        prepare_continue()
    elif opt.continue_mode == 0 and os.path.exists(opt.out_dir):
        raise Usage("Output directory " + opt.out_dir + " already exists, please change the parameter -o to another value to avoid overwriting.")

def check_opt():
    global opt
    if opt.host_mem <= 0:
        raise Usage("Please specify a positive number for -m flag.")
    elif opt.host_mem < 1:
        total_mem = detect_available_mem()
        opt.host_mem = long(total_mem * opt.host_mem)
        if total_mem <= 0:
            raise Usage("Failed to detect available memory. Please specify the value in bytes using -m flag.")
        else:
            print(str(round(total_mem/(1024**3),3)) + "Gb memory in total.", file=sys.stderr)
            print("Using: " + str(round(float(opt.host_mem)/(1024**3),3)) + "Gb.", file=sys.stderr)
    else:
        opt.host_mem = long(opt.host_mem)

    # set mode
    if opt.presets != "":
        if opt.presets == "meta":
            opt.min_count = 2
            opt.k_list = [21,41,61,81,99]
            opt.set_list_by_min_max_step = False
        elif opt.presets == "meta-sensitive":
            opt.min_count = 2
            opt.k_list = [21,31,41,51,61,71,81,91,99]
            opt.set_list_by_min_max_step = False
        elif opt.presets == "meta-large":
            opt.min_count = 2
            opt.k_list = [27,37,47,57,67,77,87]
            opt.set_list_by_min_max_step = False
        elif opt.presets == "bulk":
            opt.min_count = 3
            opt.no_mercy = True
            opt.k_list = [31,51,71,91,99]
            opt.set_list_by_min_max_step = False
        elif opt.presets == "single-cell":
            opt.min_count = 3
            opt.k_list = [21,33,55,77,99,121]
            opt.merge_similar, opt.merge_len = 0.96, 20
            opt.set_list_by_min_max_step = False
        else:
            raise Usage("Invalid preset: " + opt.presets)

    if opt.set_list_by_min_max_step:
        if opt.k_step % 2 == 1:
            raise Usage("k-step must be even number!")
        if opt.k_min > opt.k_max:
            raise Usage("Error: k_min > k_max!")

        opt.k_list = list()
        k = opt.k_min
        while k < opt.k_max:
            opt.k_list.append(k)
            k = k + opt.k_step
        opt.k_list.append(opt.k_max)

    if len(opt.k_list) == 0:
        raise Usage("k list should not be empty!")

    if opt.k_list[0] < 15 or opt.k_list[len(opt.k_list) - 1] > 127:
        raise Usage("All k's should be in range [15, 127]")

    # for k in opt.k_list:
    #     if k % 2 == 0:
    #         raise Usage("All k must be odd number!")

    for i in range(1, len(opt.k_list)):
        if opt.k_list[i] - opt.k_list[i-1] > 28:
            raise Usage("k-step/adjacent k difference must be <= 28")

    opt.k_min, opt.k_max = opt.k_list[0], opt.k_list[len(opt.k_list) - 1]

    if opt.use_gpu == 0:
        opt.gpu_mem = 0
    if opt.k_max < opt.k_min:
        raise Usage("k_min should be no larger than k_max.")
    if opt.min_count <= 0:
        raise Usage("min_count must be greater than 0.")
    elif opt.min_count == 1:
        opt.kmin_1pass = True
        opt.no_mercy = True
    if opt.prune_level < 0 or opt.prune_level > 2:
        raise Usage("prune level must be in 0-2.")
    if opt.merge_len < 0:
        raise Usage("merge_level: length must be >= 0")
    if opt.merge_similar < 0 or opt.merge_similar > 1:
        raise Usage("merge_level: similarity must be in [0, 1]")
    if opt.low_local_ratio <= 0 or opt.low_local_ratio > 0.5:
        raise Usage("low_local_ratio should be in (0, 0.5].")
    if opt.num_cpu_threads > multiprocessing.cpu_count():
        print("Maximum number of available CPU thread is %d." % multiprocessing.cpu_count(), file=sys.stderr);
        print("Number of thread is reset to the %d." % max(2, multiprocessing.cpu_count()), file=sys.stderr);
        opt.num_cpu_threads = multiprocessing.cpu_count()
    if opt.num_cpu_threads == 0:
        opt.num_cpu_threads = multiprocessing.cpu_count()
    if opt.num_cpu_threads <= 1:
        raise Usage("num_cpu_threads should be at least 2.")

    # reads
    if len(opt.pe1) != len(opt.pe2):
        raise Usage("Number of paired-end files not match!")
    for r in opt.pe1 + opt.pe2 + opt.se + opt.pe12:
        if not os.path.exists(r):
            raise Usage("Cannot find file " + r)

    if opt.input_cmd == "" and len(opt.pe1 + opt.pe2 + opt.se + opt.pe12) == 0:
        raise Usage("No input files or input command!")

def detect_available_mem():
    mem = long()
    if sys.platform.find("linux") != -1:
        try:
            mem = long(float(os.popen("free").readlines()[1].split()[1]) * 1024)
        except IndexError:
            mem = 0
    elif sys.platform.find("darwin") != -1:
        try:
            mem = long(float(os.popen("sysctl hw.memsize").readlines()[0].split()[1]))
        except IndexError:
            mem = 0
    else:
        mem = 0
    return mem

def write_opt(argv):
    with open(opt_file_name(), "w") as f:
        print("\n".join(argv), file=f)
    f.close()

def prepare_continue():
    global opt # out_dir is already set
    if not os.path.exists(opt_file_name()):
        print("Cannot find " + opt.out_dir + "opts.txt", file=sys.stderr)
        print("Please check whether the output directory is correctly set by \"-o\"", file=sys.stderr)
        print("Now switching to normal mode.", file=sys.stderr)
        return

    print("Continue mode activated. Ignore all options other than -o/--out-dir.", file=sys.stderr)

    with open(opt_file_name(), "r") as f:
        argv = []
        for line in f:
            argv.append(line.strip())
        print("Continue with options: " + " ".join(argv), file=sys.stderr)
        t_dir = opt.out_dir
        opt = Options()
        opt.out_dir = t_dir
        opt.continue_mode = True # avoid dead loop
        parse_opt(argv)
    f.close()

    opt.last_cp = -1
    if os.path.exists(opt.temp_dir + "cp.txt"):
        with open(opt.temp_dir + "cp.txt", "r") as cpf:
            for line in cpf:
                a = line.strip().split()
                if len(a) == 2 and a[1] == "done":
                    opt.last_cp = int(a[0])
        cpf.close()
    print("Continue from check point " + str(opt.last_cp), file=sys.stderr)

def check_bin():
    for subprogram in ["megahit_gt"]:
        if not os.path.exists(opt.bin_dir + subprogram):
            raise Usage("Cannot find sub-program \"" + subprogram + "\", please recompile.")

def get_version():
    global megahit_version_str
    megahit_version_str = "MEGAHIT-GT " + \
                          subprocess.Popen([opt.bin_dir + "megahit_gt", "dumpversion"],
                                           stdout=subprocess.PIPE).communicate()[0].rstrip().decode('utf-8')

def graph_prefix(kmer_k):
    if not os.path.exists(opt.temp_dir + "k" + str(kmer_k)):
        os.mkdir(opt.temp_dir + "k" + str(kmer_k))
    return opt.temp_dir + "k" + str(kmer_k) + "/" + str(kmer_k)

def delect_file_if_exist(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)

def delete_tmp_after_build(kmer_k):
    for i in range(0, opt.num_cpu_threads):
        delect_file_if_exist(graph_prefix(kmer_k) + ".edges." + str(i))
    for i in range(0, 64):
        delect_file_if_exist(graph_prefix(kmer_k) + ".mercy_cand." + str(i))
    for i in range(0, opt.num_cpu_threads - 1):
        delect_file_if_exist(graph_prefix(kmer_k) + ".mercy." + str(i))
    delect_file_if_exist(graph_prefix(kmer_k) + ".cand")

def make_out_dir():
    if os.path.exists(opt.out_dir):
        pass
    else:
        os.mkdir(opt.out_dir)

    if os.path.exists(opt.temp_dir):
        pass
    else:
        os.mkdir(opt.temp_dir)

    if os.path.exists(opt.contig_dir):
        pass
    else:
        os.mkdir(opt.contig_dir)

def write_cp():
    global cp
    cpf = open(opt.temp_dir + "cp.txt", "a")
    print(str(cp) + "\t" + "done", file=cpf);
    cp = cp + 1
    cpf.close()

def inpipe_cmd(file_name):
    if file_name.endswith('.gz'):
        return 'gzip -cd ' + file_name
    elif file_name.endswith('.bz2'):
        return 'bzip2 -cd ' + file_name
    else:
        return "cat " + file_name

def write_lib():
    global opt
    opt.lib = opt.temp_dir + "reads.lib"
    lib = open(opt.lib, "w")
    for i in range(0, len(opt.pe12)):
        print(opt.pe12[i], file=lib)

        if inpipe_cmd(opt.pe12[i]) != "":
            print("interleaved " + opt.temp_dir + "inpipe.pe12." + str(i), file=lib)
        else:
            print("interleaved " + opt.pe12[i], file=lib)

    for i in range(0, len(opt.pe1)):

        if inpipe_cmd(opt.pe1[i]) != "":
            f1 = opt.temp_dir + "inpipe.pe1." + str(i)
        else:
            f1 = opt.pe1[i]

        if inpipe_cmd(opt.pe2[i]) != "":
            f2 = opt.temp_dir + "inpipe.pe2." + str(i)
        else:
            f2 = opt.pe2[i]

        print(','.join([opt.pe1[i], opt.pe2[i]]), file=lib)
        print("pe " + f1 + " " + f2, file=lib)

    for i in range(0, len(opt.se)):
        print(opt.se[i], file=lib)

        if inpipe_cmd(opt.se[i]) != "":
            print("se " + opt.temp_dir + "inpipe.se." + str(i), file=lib)
        else:
            print("se " + opt.se[i], file=lib)

    if opt.input_cmd != "":
        print('\"' + opt.input_cmd + '\"', file=lib)
        print("se " + "-", file=lib)

    lib.close()

def build_lib():
    global cp
    if (not opt.continue_mode) or (cp > opt.last_cp):
        build_lib_cmd = [opt.bin_dir + "megahit_gt", "buildlib",
                         opt.lib,
                         opt.lib]

        fifos = list()
        pipes = list()
        try:
            # create inpipe

            for i in range(0, len(opt.pe12)):
                if inpipe_cmd(opt.pe12[i]) != "":
                    delect_file_if_exist(opt.temp_dir + "inpipe.pe12." + str(i))
                    os.mkfifo(opt.temp_dir + "inpipe.pe12." + str(i))
                    fifos.append(opt.temp_dir + "inpipe.pe12." + str(i))

            for i in range(0, len(opt.pe1)):
                if inpipe_cmd(opt.pe1[i]) != "":
                    delect_file_if_exist(opt.temp_dir + "inpipe.pe1." + str(i))
                    os.mkfifo(opt.temp_dir + "inpipe.pe1." + str(i))
                    fifos.append(opt.temp_dir + "inpipe.pe1." + str(i))
                
                if inpipe_cmd(opt.pe2[i]) != "":
                    delect_file_if_exist(opt.temp_dir + "inpipe.pe2." + str(i))
                    os.mkfifo(opt.temp_dir + "inpipe.pe2." + str(i))
                    fifos.append(opt.temp_dir + "inpipe.pe2." + str(i))

            for i in range(0, len(opt.se)):
                if inpipe_cmd(opt.se[i]) != "":
                    delect_file_if_exist(opt.temp_dir + "inpipe.se." + str(i))
                    os.mkfifo(opt.temp_dir + "inpipe.se." + str(i))
                    fifos.append(opt.temp_dir + "inpipe.se." + str(i))

            logging.info("--- [%s] Converting reads to binaries ---" % datetime.now().strftime("%c"))
            logging.debug("%s" % (" ").join(build_lib_cmd))

            if opt.input_cmd != "":
                logging.debug("input cmd: " + opt.input_cmd)
                input_thread = subprocess.Popen(opt.input_cmd, shell = True, stdout = subprocess.PIPE)
                p = subprocess.Popen(build_lib_cmd, stdin = input_thread.stdout, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            else:
                p = subprocess.Popen(build_lib_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # output to inpipe

            for i in range(0, len(opt.pe12)):
                if inpipe_cmd(opt.pe12[i]) != "":
                    ip_thread12 = subprocess.Popen(inpipe_cmd(opt.pe12[i]) + " > " + opt.temp_dir + "inpipe.pe12." + str(i), shell = True, preexec_fn = os.setsid)
                    pipes.append(ip_thread12)

            for i in range(0, len(opt.pe1)):
                if inpipe_cmd(opt.pe1[i]) != "":
                    ip_thread1 = subprocess.Popen(inpipe_cmd(opt.pe1[i]) + " > " + opt.temp_dir + "inpipe.pe1." + str(i), shell = True, preexec_fn = os.setsid)
                    pipes.append(ip_thread1)
                
                if inpipe_cmd(opt.pe2[i]) != "":
                    ip_thread2 = subprocess.Popen(inpipe_cmd(opt.pe2[i]) + " > " + opt.temp_dir + "inpipe.pe2." + str(i), shell = True, preexec_fn = os.setsid)
                    pipes.append(ip_thread2)

            for i in range(0, len(opt.se)):
                if inpipe_cmd(opt.se[i]) != "":
                    ip_thread_se = subprocess.Popen(inpipe_cmd(opt.se[i]) + " > " + opt.temp_dir + "inpipe.se." + str(i), shell = True, preexec_fn = os.setsid)
                    pipes.append(ip_thread_se)
            
            while True:
                line = p.stderr.readline().rstrip()
                if not line:
                    break;
                logging.info(line)

            ret_code = p.wait()

            if ret_code != 0:
                logging.error("Error occurs when running \"megahit_asm_core buildlib\"; please refer to %s for detail" % log_file_name())
                logging.error("[Exit code %d]" % ret_code)
                exit(ret_code)

            while len(pipes) > 0:
                pp = pipes.pop()
                pp_ret = pp.wait()
                if pp_ret != 0:
                    logging.error("Error occurs when reading inputs")
                    exit(pp_ret)

        except OSError as o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                logging.error("Error: sub-program megahit_asm_core not found, please recompile MEGAHIT")
            exit(1)
        except KeyboardInterrupt:
            p.terminate()
            exit(1)

        finally:
            for p in pipes:
                os.killpg(p.pid, signal.SIGTERM)
            for f in fifos:
                delect_file_if_exist(f)

    write_cp()

def build_graph(list_of_k):
    global cp
    phase1_out_threads = max(1, int(opt.num_cpu_threads / 3))
    if (not opt.continue_mode) or (cp > opt.last_cp):
        count_opt = ["-k", str(opt.k_current),
                     "-m", str(opt.min_count),
                     "--host_mem", str(opt.host_mem),
                     "--mem_flag", str(opt.mem_flag),
                     "--gpu_mem", str(opt.gpu_mem),
                     "--output_prefix", graph_prefix(opt.k_current),
                     "--num_cpu_threads", str(opt.num_cpu_threads),
                     "--num_output_threads", str(phase1_out_threads),
                     "--read_lib_file", opt.lib]

        cmd = [opt.bin_dir + "megahit_gt", "buildgraph"] + count_opt
        if not opt.no_mercy:
            cmd.append("--need_mercy")

        if list_of_k and opt.k_current!=opt.k_min:
            previous_k = opt.k_list[opt.k_list.index(opt.k_current) - 1]
            # a set of contigs is used to assist the next graph, how?
            cmd += ["--assist_seq", opt.combined_contigs_file]

        try:
            if opt.kmin_1pass:
                logging.info("--- [%s] Extracting solid (k+1)-mers and building sdbg for k = %d ---" % (datetime.now().strftime("%c"), opt.k_current))
            else:
                logging.info("--- [%s] Extracting solid (k+1)-mers for k = %d ---" % (datetime.now().strftime("%c"), opt.k_current))

            logging.debug("cmd: %s" % (" ").join(cmd))
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            while True:
                line = p.stderr.readline().rstrip()
                if not line:
                    break;
                logging.debug(line)

            ret_code = p.wait()

            if ret_code != 0:
                logging.error("Error occurs when running \"sdbg_builder count/read2sdbg\", please refer to %s for detail" % log_file_name())
                logging.error("[Exit code %d] " % ret_code)
                exit(ret_code)

        except OSError as o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                logging.error("Error: sub-program sdbg_builder not found, please recompile MEGAHIT-GT")
            exit(1)
        except KeyboardInterrupt:
            p.terminate()
            exit(1)

    write_cp()

def parse_gene_list():
    global cp
    with open(opt.gene_list, "r") as f:
        for line in f:
            words = line.split()
            # print (words)
            opt.gene_info[words[0]] = [words[1], words[2], words[3]]
            # print (opt.gene_info)
    f.close()

def find_seed(gene):
    global cp
    parameter = [opt.gene_info[gene][2], str(opt.se[0]), str(opt.k_current + 1), str(opt.num_cpu_threads)]
    cmd = [opt.bin_dir + "megahit_gt", "findstart"] + parameter

    try:
        logging.info("--- [%s] Finding starting kmers for %s k = %d ---" % (datetime.now().strftime("%c"), gene, opt.k_current))
        logging.debug("cmd: %s" % (" ").join(cmd))
        seed_file = opt.out_dir + str(opt.k_current) + "_" + gene + "_starting_kmers_unsorted.txt"
        final_seed_file = opt.out_dir + str(opt.k_current) + "_" + gene + "_starting_kmers.txt"
        with open(seed_file, "w") as starting_kmers:
            p = subprocess.Popen(cmd, stdout=starting_kmers, stderr=subprocess.PIPE)
        ret_code = p.wait()
        starting_kmers.close()

        if ret_code == 0:
            with open(final_seed_file, "w") as final_starting_kmers:
                sort_cmd = "sort -k4 " + seed_file + " | uniq | sort -nk8"
                p = subprocess.Popen(sort_cmd, shell = True, stdout=final_starting_kmers, stderr=subprocess.PIPE)

    except OSError as o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            logging.error("Error: sub-program megahit_gt not found, please recompile MEGAHIT-GT")
        exit(1)
    except KeyboardInterrupt:
        p.terminate()
        exit(1)
    write_cp()

def search_contigs():
    global cp
    parameter = [graph_prefix(opt.k_current), opt.gene_list, opt.out_dir + str(opt.k_current), opt.out_dir + str(opt.k_current), str(opt.num_cpu_threads)]
    cmd = [opt.bin_dir + "megahit_gt", "search"] + parameter

    try:
        logging.info("--- [%s] Searching contigs for k = %d ---" % (datetime.now().strftime("%c"), opt.k_current))
        logging.debug("cmd: %s" % (" ").join(cmd))
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE)

        while True:
            line = p.stderr.readline().rstrip()
            if not line:
                break;
            logging.debug(line)
        p.wait()

        cmd2 = "cat " + opt.out_dir + str(opt.k_current) + "_raw_contigs*"
        logging.debug("cmd: %s" % (" ").join(cmd2))
        opt.combined_contigs_file = opt.out_dir + str(opt.k_current) + "_combined_contigs" + ".fasta"
        with open(opt.combined_contigs_file, "w") as raw_contigs:
            p = subprocess.Popen(cmd2, shell = True, stdout=raw_contigs, stderr=subprocess.PIPE)
        ret_code = p.wait()
        raw_contigs.close()

        if ret_code == 0:
            cmd3 = "cat " + opt.combined_contigs_file + " | " +  opt.bin_dir + "megahit_gt readstat" + "| head -n 2 | cut -d ' ' -f  3 | paste -d ' ' -s"
            logging.debug(cmd3)
            stat_file = opt.combined_contigs_file + ".info"
            with open(stat_file, "w") as raw_contigs_info:
                p = subprocess.Popen(cmd3, shell = True, stdout=raw_contigs_info, stderr=subprocess.PIPE)
            p.wait()
    except OSError as o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            logging.error("Error: sub-program megahit_gt not found, please recompile MEGAHIT-GT")
        exit(1)

    except KeyboardInterrupt:
        p.terminate()
        exit(1)
    write_cp()

def filter_contigs(input_file, output_file):
    global cp
    parameter = [str(opt.filter_len)]
    cmd = [opt.bin_dir + "megahit_gt", "filterbylen"] + parameter

    try:
        logging.info("--- [%s] Filtering contigs with min_len = %d ---" % (datetime.now().strftime("%c"), opt.filter_len))
        logging.debug("cmd: %s" % (" ").join(cmd))
        # nucl_file = opt.out_dir + opt.filtered_nucl_file
        with open(output_file, "w") as filtered_nucl_contigs:
            with open(input_file, "r") as raw_contigs:
                p = subprocess.Popen(cmd, stdin=raw_contigs, stdout=filtered_nucl_contigs, stderr=subprocess.PIPE)
        p.wait()
    except OSError as o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            logging.error("Error: sub-program megahit_toolkit not found, please recompile MEGAHIT")
        exit(1)
    except KeyboardInterrupt:
        p.terminate()
        exit(1)
    write_cp()

def translate_to_aa(input_file, output_file):
    global cp
    cmd = [opt.bin_dir + "megahit_gt", "translate"] + input_file

    try:
        logging.info("--- [%s] Translating nucl contigs to aa contigs ---" % (datetime.now().strftime("%c")))
        logging.debug("cmd: %s" % (" ").join(cmd))
        with open(output_file, "w") as filtered_prot_contigs:
            p = subprocess.Popen(cmd, stdout=filtered_prot_contigs, stderr=subprocess.PIPE)
        p.wait()
    except OSError as o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            logging.error("Error: sub-program translate not found, please recompile MEGAHIT-GT")
        exit(1)
    except KeyboardInterrupt:
        p.terminate()
        exit(1)
    write_cp()

def post_processing():
    global cp
    try:
        logging.info("--- [%s] Preparing for post-processing ---" % (datetime.now().strftime("%c")))
        post_proc_directory = opt.out_dir + "for_post_proc/"
        os.makedirs(post_proc_directory)
        for key in opt.gene_info:
            os.makedirs(post_proc_directory + key)
            # shutil.copyfile(opt.out_dir+str(opt.k_current)+"_raw_contigs_"+key+".fasta", post_proc_directory+key+"/nucl_merged.fasta")
            filter_contigs(opt.out_dir+str(opt.k_current)+"_raw_contigs_"+key+".fasta", post_proc_directory+key+"/nucl_merged.fasta")
            translate_to_aa([post_proc_directory+key+"/nucl_merged.fasta"], post_proc_directory+key+"/prot_merged.fasta")
    except OSError as o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            logging.error("Error: XXXX, please recompile MEGAHIT-GT")
        exit(1)
    except KeyboardInterrupt:
        p.terminate()
        exit(1)
    write_cp()

def main(argv = None):
    if argv is None:
        argv = sys.argv

    try:
        start_time = time.time()

        check_bin()
        get_version()
        parse_opt(argv[1:])
        check_opt()
        make_out_dir()

        logging.basicConfig(level = logging.NOTSET,
                format = '%(message)s',
                filename = log_file_name(),
                filemode = 'a')

        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        if opt.verbose:
            console.setLevel(logging.NOTSET)

        formatter = logging.Formatter('%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

        logging.info(megahit_version_str)
        logging.info("--- [%s] Start assembly. Number of CPU threads %d ---" % (datetime.now().strftime("%c"), opt.num_cpu_threads))
        logging.info("--- [%s] k list: %s ---" % (datetime.now().strftime("%c"), ','.join(map(str, opt.k_list))))

        if not opt.continue_mode:
            write_opt(argv[1:]) # for --continue

        write_lib()
        build_lib()
        parse_gene_list()

        if len(opt.k_list) == 1:
            opt.k_current = opt.k_list[0]
            build_graph(False)
            for key in opt.gene_info:
                find_seed(key)
            search_contigs()
            post_processing()
        elif len(opt.k_list) > 1:
            for k in opt.k_list:
                opt.k_current = k
                build_graph(True)
                # for each gene
                for key in opt.gene_info:
                    find_seed(key)
                search_contigs()
            post_processing()


        logging.info("--- [%s] ALL DONE. Time elapsed: %f seconds ---" % (datetime.now().strftime("%c"), time.time() - start_time))

    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg), file=sys.stderr)
        return 2

if __name__ == "__main__":
    sys.exit(main())
