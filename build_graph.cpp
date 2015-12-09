/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include <stdio.h>
#include <omp.h>

#include <iostream>
#include <string>
#include <stdexcept>

#include "cx1_read2sdbg.h"
#include "options_description.h"
#include "utils.h"
#include "definitions.h"

int build_graph(int argc, char **argv) {
    AutoMaxRssRecorder recorder;

    // parse option the same as kmer_count
    OptionsDescription desc;
    read2sdbg_opt_t opt;

    desc.AddOption("kmer_k", "k", opt.kmer_k, "kmer size");
    desc.AddOption("min_kmer_frequency", "m", opt.kmer_freq_threshold, "min frequency to output an edge");
    desc.AddOption("host_mem", "", opt.host_mem, "Max memory to be used. 90% of the free memory is recommended.");
    desc.AddOption("gpu_mem", "", opt.gpu_mem, "gpu memory to be used. 0 for auto detect.");
    desc.AddOption("num_cpu_threads", "", opt.num_cpu_threads, "number of CPU threads. At least 2.");
    desc.AddOption("num_output_threads", "", opt.num_output_threads, "number of threads for output. Must be less than num_cpu_threads");
    desc.AddOption("read_lib_file", "", opt.read_lib_file, "input fast[aq] file, can be gzip'ed. \"-\" for stdin.");
    desc.AddOption("assist_seq", "", opt.assist_seq_file, "input assisting fast[aq] file (FILE_NAME.info should exist), can be gzip'ed.");
    desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
    desc.AddOption("mem_flag", "", opt.mem_flag, "memory options. 0: minimize memory usage; 1: automatically use moderate memory; other: use all available mem specified by '--host_mem'");
    desc.AddOption("need_mercy", "", opt.need_mercy, "to add mercy edges.");

    try {
        desc.Parse(argc, argv);

        if (opt.read_lib_file == "") {
            throw std::logic_error("No input file!");
        }

        if (opt.num_cpu_threads == 0) {
            opt.num_cpu_threads = omp_get_max_threads();
        }

        if (opt.num_output_threads == 0) {
            opt.num_output_threads = std::max(1, opt.num_cpu_threads / 3);
        }

        if (opt.host_mem == 0) {
            throw std::logic_error("Please specify the host memory!");
        }

        if (opt.num_cpu_threads == 1) {
            throw std::logic_error("Number of CPU threads should be at least 2!");
        }

        if (opt.num_output_threads >= opt.num_cpu_threads) {
            throw std::logic_error("Number of output threads must be less than number of CPU threads!");
        }
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: sdbg_builder read2sdbg --read_lib_file fastx_file -o out" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }

    cx1_read2sdbg::read2sdbg_global_t globals;
    globals.kmer_k = opt.kmer_k;
    globals.kmer_freq_threshold = opt.kmer_freq_threshold;
    globals.host_mem = opt.host_mem;
    globals.gpu_mem = opt.gpu_mem;
    globals.num_cpu_threads = opt.num_cpu_threads;
    globals.num_output_threads = opt.num_output_threads;
    globals.read_lib_file = opt.read_lib_file;
    globals.assist_seq_file = opt.assist_seq_file;
    globals.output_prefix = opt.output_prefix;
    globals.mem_flag = opt.mem_flag;
    globals.need_mercy = opt.need_mercy;
    globals.cx1.g_ = &globals;

    // stage1
    if (opt.kmer_freq_threshold > 1) {
        globals.cx1.encode_lv1_diff_base_func_ = cx1_read2sdbg::s1::s1_encode_lv1_diff_base;
        globals.cx1.prepare_func_ = cx1_read2sdbg::s1::s1_read_input_prepare;
        globals.cx1.lv0_calc_bucket_size_func_ = cx1_read2sdbg::s1::s1_lv0_calc_bucket_size;
        globals.cx1.init_global_and_set_cx1_func_ = cx1_read2sdbg::s1::s1_init_global_and_set_cx1;
        globals.cx1.lv1_fill_offset_func_ = cx1_read2sdbg::s1::s1_lv1_fill_offset;
        globals.cx1.lv1_sort_and_proc = cx1_read2sdbg::s1::s1_lv1_direct_sort_and_count;
        globals.cx1.lv2_extract_substr_func_ = cx1_read2sdbg::s1::s1_lv2_extract_substr;
        globals.cx1.lv2_sort_func_ = cx1_read2sdbg::s1::s1_lv2_sort;
        globals.cx1.lv2_pre_output_partition_func_ = cx1_read2sdbg::s1::s1_lv2_pre_output_partition;
        globals.cx1.lv2_output_func_ = cx1_read2sdbg::s1::s1_lv2_output;
        globals.cx1.lv2_post_output_func_ = cx1_read2sdbg::s1::s1_lv2_post_output;
        globals.cx1.post_proc_func_ = cx1_read2sdbg::s1::s1_post_proc;
        globals.cx1.run();
    }
    else {
        cx1_read2sdbg::s1::s1_read_input_prepare(globals);
    }

    // stage2
    globals.cx1.encode_lv1_diff_base_func_ = cx1_read2sdbg::s2::s2_encode_lv1_diff_base;
    globals.cx1.prepare_func_ = cx1_read2sdbg::s2::s2_read_mercy_prepare;
    globals.cx1.lv0_calc_bucket_size_func_ = cx1_read2sdbg::s2::s2_lv0_calc_bucket_size;
    globals.cx1.init_global_and_set_cx1_func_ = cx1_read2sdbg::s2::s2_init_global_and_set_cx1;
    globals.cx1.lv1_fill_offset_func_ = cx1_read2sdbg::s2::s2_lv1_fill_offset;
    globals.cx1.lv1_sort_and_proc = cx1_read2sdbg::s2::s2_lv1_direct_sort_and_proc;
    globals.cx1.lv2_extract_substr_func_ = cx1_read2sdbg::s2::s2_lv2_extract_substr;
    globals.cx1.lv2_sort_func_ = cx1_read2sdbg::s2::s2_lv2_sort;
    globals.cx1.lv2_pre_output_partition_func_ = cx1_read2sdbg::s2::s2_lv2_pre_output_partition;
    globals.cx1.lv2_output_func_ = cx1_read2sdbg::s2::s2_lv2_output;
    globals.cx1.lv2_post_output_func_ = cx1_read2sdbg::s2::s2_lv2_post_output;
    globals.cx1.post_proc_func_ = cx1_read2sdbg::s2::s2_post_proc;
    globals.cx1.run();

    return 0;
}
