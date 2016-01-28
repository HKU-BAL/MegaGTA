#!/usr/bin/env python

#-------------------------------------------------------------------------
# MEGAHIT
# Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------

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

megahit_version_str = ""
usage_message = '''
Copyright (c) The University of Hong Kong

Usage:
  megahit_gt.py [options] {-1 <pe1> -2 <pe2> | --12 <pe12> | -r <se>} -g gene_list.txt [-o <out_dir>]

  Input options that can be specified for multiple times (supporting plain text and gz/bz2 extensions)
    -1                       <pe1>          comma-separated list of fasta/q paired-end #1 files, paired with files in <pe2>
    -2                       <pe2>          comma-separated list of fasta/q paired-end #2 files, paired with files in <pe1>
    --12                     <pe12>         comma-separated list of interleaved fasta/q paired-end files
    -r/--read                <se>           comma-separated list of fasta/q single-end files
    
    -g/--gene-list           <string>       gene list

Optional Arguments:
  Basic assembly options:
    -c/--min-count           <int>          minimum multiplicity for filtering k-mers [1]
    -k/--k-list              <int,int,..>   comma-separated list of kmer size (in range 15-63) [29,35,44]
    -p/--prune-len           <int>          prune the search if the score does not improve after <int> steps [20]
    -l/--low-cov-penalty     <float>        penalty for coverage one edges (in [0,1]) [0.5]
    --max-tip-len            <int>          max tip length [150]
    --no-mercy                              do not add mercy kmers

  Hardware options:
    -m/--memory              <float>        max memory in byte to be used in SdBG construction [0.9]
                                            (if set between 0-1, fraction of the machine's total memory)
    --mem-flag               <int>          SdBG builder memory mode [1]
                                            0: minimum; 1: moderate; others: use all memory specified by '-m/--memory'.
    --use-gpu                               use GPU
    --gpu-mem                <float>        GPU memory in byte to be used. Default: auto detect to use up all free GPU memory [0]
    -t/--num-cpu-threads     <int>          number of CPU threads, at least 2. Default: auto detect to use all CPU threads [auto]

  Output options:
    -o/--out-dir             <string>       output directory [./megahit_gt_out]
    --min-contig-len         <int>          minimum length of contigs to output [450]
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
        self.out_dir = "./megahit_gt_out/"
        self.min_contig_len = 450
        self.max_tip_len = 150
        self.prune_len = 20
        self.low_cov_penalty = 0.5
        self.k_list = [29,35,44]
        self.min_count = 1
        self.bin_dir = sys.path[0] + "/"
        self.no_mercy = False
        self.num_cpu_threads = 0
        self.temp_dir = self.out_dir + "tmp/"
        self.keep_tmp_files = False
        self.use_gpu = False
        self.mem_flag = 1
        self.continue_mode = False;
        self.last_cp = -1;
        self.kmin_1pass = False
        self.verbose = False
        self.pe1 = []
        self.pe2 = []
        self.pe12 = []
        self.aligned_ref = []
        self.se = []
        self.input_cmd = ""
        self.gene_info = {}
        self.gene_list = ""

opt = Options()
cp = 0

def log_file_name(): return opt.out_dir + "log"

def opt_file_name():
    return opt.out_dir + "opts.txt"

def mkdir(d):
    if os.path.exists(d):
        pass
    else:
        os.mkdir(d)

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
        opts, args = getopt.getopt(argv, "hm:o:r:t:v1:2:l:k:l:p:c:g:", 
                ["help",
                    "read=",
                    "12=",
                    "input-cmd=",
                    "memory=",
                    "out-dir=",
                    "min-contig-len=",
                    "max-tip-len=",
                    "--low-cov-penalty=",
                    "--prune-len=",
                    "use-gpu",
                    "num-cpu-threads=",
                    "gpu-mem=",
                    "kmin-1pass",
                    "k-list=",
                    "num-cpu-threads=",
                    "min-count=",
                    "no-mercy",
                    "keep-tmp-files",
                    "mem-flag=",
                    "continue",
                    "version",
                    "verbose",
                    "gene-list="])
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
        elif option in ("--k-list", "-k"):
            opt.k_list = list(map(int, value.split(",")))
            opt.k_list.sort()
        elif option in ("--min-count", "-c"):
            opt.min_count = int(value)
        elif option == "--max-tip-len":
            opt.max_tip_len = int(value)
        elif option == "--no-mercy":
            opt.no_mercy = True
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
        elif option in ("-r", "--read"):
            opt.se += value.split(",")
        elif option == "-1":
            opt.pe1 += value.split(",")
        elif option == "-2":
            opt.pe2 += value.split(",")
        elif option == "--12":
            opt.pe12 += value.split(",")
        elif option in ("--gene-list", "-g"):
            opt.gene_list = value
        elif option in ("--prune-len", "-p"):
            opt.prune_len = int(value)
        elif option in ("--low-cov-penalty", "-l"):
            opt.low_cov_penalty = float(value)

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

    if len(opt.k_list) == 0:
        raise Usage("k list should not be empty!")

    if opt.k_list[0] < 15 or opt.k_list[len(opt.k_list) - 1] > 127:
        raise Usage("All k's should be in range [15, 127]")

    if opt.use_gpu == 0:
        opt.gpu_mem = 0
    if opt.min_count <= 0:
        raise Usage("min_count must be greater than 0.")
    elif opt.min_count == 1:
        opt.kmin_1pass = True
        opt.no_mercy = True
    if opt.num_cpu_threads > multiprocessing.cpu_count():
        print("Maximum number of available CPU thread is %d." % multiprocessing.cpu_count(), file=sys.stderr);
        print("Number of thread is reset to the %d." % max(2, multiprocessing.cpu_count()), file=sys.stderr);
        opt.num_cpu_threads = multiprocessing.cpu_count()
    if opt.num_cpu_threads == 0:
        opt.num_cpu_threads = multiprocessing.cpu_count()
    if opt.gene_list == "":
        raise Usage("--gene-list could not be empty")
    if opt.prune_len <= 0:
        raise Usage("prune length should be >= 1")
    if opt.low_cov_penalty < 0 or opt.low_cov_penalty > 1:
        raise Usage("low coverage penalty should be between [0, 1]")

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
    if not os.path.exists(opt.out_dir + "k" + str(kmer_k)):
        os.mkdir(opt.out_dir + "k" + str(kmer_k))
    return opt.out_dir + "k" + str(kmer_k) + "/" + str(kmer_k)

def contig_file(kmer_k):
    return graph_prefix(kmer_k) + ".contigs.fa"

def delect_file_if_exist(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)

def make_out_dir():
    mkdir(opt.out_dir)
    mkdir(opt.temp_dir)

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

def build_graph(k, assist_seq):
    global cp
    phase1_out_threads = max(1, int(opt.num_cpu_threads / 3))
    if (not opt.continue_mode) or (cp > opt.last_cp):
        count_opt = ["-k", str(k),
                     "-m", str(opt.min_count),
                     "--host_mem", str(opt.host_mem),
                     "--mem_flag", str(opt.mem_flag),
                     "--gpu_mem", str(opt.gpu_mem),
                     "--output_prefix", graph_prefix(k),
                     "--num_cpu_threads", str(opt.num_cpu_threads),
                     "--num_output_threads", str(phase1_out_threads),
                     "--read_lib_file", opt.lib]

        cmd = [opt.bin_dir + "megahit_gt", "buildgraph"] + count_opt
        if not opt.no_mercy:
            cmd.append("--need_mercy")

        if assist_seq != "":
            cmd += ["--assist_seq", assist_seq]

        try:
            logging.info("--- [%s] Building sdbg for k = %d ---" % (datetime.now().strftime("%c"), k))

            logging.debug("cmd: %s" % (" ").join(cmd))
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            while True:
                line = p.stderr.readline().rstrip()
                if not line:
                    break;
                logging.debug(line)

            ret_code = p.wait()

            if ret_code != 0:
                logging.error("Error occurs when running \"megahit_gt count/read2sdbg\", please refer to %s for detail" % log_file_name())
                logging.error("[Exit code %d] " % ret_code)
                exit(ret_code)

        except OSError as o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                logging.error("Error: sub-program megahit_gt not found, please recompile MEGAHIT-GT")
            exit(1)
        except KeyboardInterrupt:
            p.terminate()
            exit(1)

    write_cp()

def assemble(k):
    global cp
    if (not opt.continue_mode) or (cp > opt.last_cp):
        min_standalone = 400 # TODO HARDCODE

        assembly_cmd = [opt.bin_dir + "megahit_gt", "denovo",
                        "-s", graph_prefix(k),
                        "-o", graph_prefix(k),
                        "-t", str(opt.num_cpu_threads),
                        "--min_standalone", str(min_standalone),
                        "--max_tip_len", str(opt.max_tip_len)]

        index_k = opt.k_list.index(k);
        assembly_cmd += ["--min_contig", str(opt.k_list[index_k + 1] + 1)]
            
        try:
            logging.info("--- [%s] De novo assembling contigs from SdBG for k = %d ---" % (datetime.now().strftime("%c"), k))
            logging.debug("cmd: %s" % (" ").join(assembly_cmd))

            p = subprocess.Popen(assembly_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            
            while True:
                line = p.stderr.readline().rstrip()
                if not line:
                    break;
                logging.debug(line)

            ret_code = p.wait()

            if ret_code != 0:
                logging.error("Error occurs when assembling contigs for k = %d, please refer to %s for detail" % (k, log_file_name()))
                logging.error("[Exit code %d]" % ret_code)
                exit(ret_code)
            
        except OSError as o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                logging.error("Error: sub-program megahit_assemble not found, please recompile MEGAHIT")
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
            opt.gene_info[words[0]] = [words[1], words[2], words[3]]
    f.close()

def find_seed(k, gene):
    global cp
    if (not opt.continue_mode) or (cp > opt.last_cp):
        parameter = [opt.gene_info[gene][2], str(opt.lib + ".bin"), str(k + 1), str(opt.num_cpu_threads)]
        # index_k = opt.k_list.index(k)
        # if index_k > 0:
        #     parameter += [contig_file(opt.k_list[index_k - 1])]
        cmd = [opt.bin_dir + "megahit_gt", "findstart"] + parameter

        try:
            logging.info("--- [%s] Finding starting kmers for %s k = %d ---" % (datetime.now().strftime("%c"), gene, k))
            logging.debug("cmd: %s" % (" ").join(cmd))
            final_seed_file = graph_prefix(k) + "_" + gene + "_starting_kmers.txt"
            with open(final_seed_file, "w") as starting_kmers:
                p = subprocess.Popen(cmd, stdout = starting_kmers, stderr = subprocess.PIPE)

                while True:
                    line = p.stderr.readline().rstrip()
                    if not line:
                        break;
                    logging.debug(line)

            ret_code = p.wait()
            starting_kmers.close()

            if ret_code != 0:
                logging.error("Error occurs when finding seeds for k = %d, please refer to %s for detail" % (k, log_file_name()))
                logging.error("[Exit code %d]" % ret_code)
                exit(ret_code)

        except OSError as o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                logging.error("Error: sub-program megahit_gt not found, please recompile MEGAHIT-GT")
            exit(1)
        except KeyboardInterrupt:
            p.terminate()
            exit(1)
    write_cp()

def search_contigs(k):
    global cp
    if (not opt.continue_mode) or (cp > opt.last_cp):
        parameter = [graph_prefix(k), opt.gene_list, graph_prefix(k), graph_prefix(k),
                     str(opt.prune_len), str(opt.low_cov_penalty), str(min(6, opt.num_cpu_threads))]
        cmd = [opt.bin_dir + "megahit_gt", "search"] + parameter

        try:
            logging.info("--- [%s] Searching contigs for k = %d ---" % (datetime.now().strftime("%c"), k))
            logging.debug("cmd: %s" % (" ").join(cmd))
            p = subprocess.Popen(cmd, stderr=subprocess.PIPE)

            while True:
                line = p.stderr.readline().rstrip()
                if not line:
                    break;
                logging.debug(line)

            ret_code = p.wait()

            if ret_code != 0:
                logging.error("Error occurs when searching contigs for k = %d, please refer to %s for detail" % (k, log_file_name()))
                logging.error("[Exit code %d]" % ret_code)
                exit(ret_code)

            # do this before combining
            post_proc_directory = opt.out_dir + "contigs/"
            mkdir(post_proc_directory)
            for gene_name in opt.gene_info:
                mkdir(post_proc_directory + gene_name)
                filter_contigs(graph_prefix(k) + "_raw_contigs_" + gene_name + ".fasta", post_proc_directory + gene_name + "/nucl_merged.fasta")
                translate_to_aa(post_proc_directory + gene_name + "/nucl_merged.fasta", post_proc_directory + gene_name + "/prot_merged.fasta")

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
    if (not opt.continue_mode) or (cp > opt.last_cp):
        parameter = [str(opt.min_contig_len)]
        cmd = [opt.bin_dir + "megahit_gt", "filterbylen"] + parameter

        try:
            logging.info("--- [%s] Filtering contigs with minimum length = %d %s->%s ---" % (datetime.now().strftime("%c"), opt.min_contig_len, input_file, output_file))
            logging.debug("cmd: %s" % (" ").join(cmd))
            with open(output_file, "w") as filtered_nucl_contigs:
                with open(input_file, "r") as raw_contigs:
                    p = subprocess.Popen(cmd, stdin = raw_contigs, stdout = filtered_nucl_contigs, stderr = subprocess.PIPE)
            ret_code = p.wait()

            if ret_code != 0:
                logging.error("Error occurs when filtering contigs for k = %d, please refer to %s for detail" % (k, log_file_name()))
                logging.error("[Exit code %d]" % ret_code)
                exit(ret_code)

        except OSError as o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                logging.error("Error: sub-program megahit_toolkit not found, please recompile MEGAHIT-GT")
            exit(1)

        except KeyboardInterrupt:
            p.terminate()
            exit(1)
    write_cp()

def translate_to_aa(input_file, output_file):
    global cp
    if (not opt.continue_mode) or (cp > opt.last_cp):
        cmd = [opt.bin_dir + "megahit_gt", "translate"] + [input_file]

        try:
            logging.info("--- [%s] Translating nucl contigs to aa contigs %s->%s ---" % (datetime.now().strftime("%c"), input_file, output_file))
            logging.debug("cmd: %s" % (" ").join(cmd))
            with open(output_file, "w") as filtered_prot_contigs:
                p = subprocess.Popen(cmd, stdout = filtered_prot_contigs, stderr = subprocess.PIPE)
            ret_code = p.wait()

            if ret_code != 0:
                logging.error("Error occurs when translating contigs for k = %d, please refer to %s for detail" % (k, log_file_name()))
                logging.error("[Exit code %d]" % ret_code)
                exit(ret_code)

        except OSError as o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                logging.error("Error: sub-program translate not found, please recompile MEGAHIT-GT")
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
        logging.info("--- [%s] Start. Number of CPU threads %d ---" % (datetime.now().strftime("%c"), opt.num_cpu_threads))
        logging.info("--- [%s] k list: %s ---" % (datetime.now().strftime("%c"), ','.join(map(str, opt.k_list))))

        if not opt.continue_mode:
            write_opt(argv[1:]) # for --continue

        write_lib()
        build_lib()
        parse_gene_list()

        for i in range(len(opt.k_list)):
            k = opt.k_list[i]
            assist_seq = ""
            if i > 0:
                assist_seq = contig_file(opt.k_list[i-1])
            build_graph(k, assist_seq)

            if i != (len(opt.k_list)) - 1:
                assemble(k)
            else:
                for gene_name in opt.gene_info:
                    find_seed(k, gene_name)
                search_contigs(k)

        logging.info("--- [%s] ALL DONE. Time elapsed: %f seconds ---" % (datetime.now().strftime("%c"), time.time() - start_time))

    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg), file=sys.stderr)
        return 2

if __name__ == "__main__":
    sys.exit(main())
