/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
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

#include <omp.h>
#include <assert.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "succinct_dbg.h"
#include "assembly_algorithms.h"
#include "utils.h"
#include "options_description.h"
#include "mem_file_checker-inl.h"
#include "unitig_graph.h"
#include "histgram.h"

using std::string;

struct asm_opt_t {
    string sdbg_name;
    string output_prefix;
    int num_cpu_threads;
    int max_tip_len;
    bool no_bubble;
    int min_standalone;
    int min_contig;
    
    asm_opt_t() {
        output_prefix = "out";
        num_cpu_threads = 0;
        max_tip_len = -1;
        no_bubble = false;
        min_standalone = 400;
        min_contig = 0;
    }

    string contig_file() {
        return output_prefix + ".contigs.fa";
    }
};

static asm_opt_t opt;

void ParseAsmOption(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("sdbg_name", "s", opt.sdbg_name, "succinct de Bruijn graph name");
    desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
    desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads, "number of cpu threads");
    desc.AddOption("max_tip_len", "", opt.max_tip_len, "max length for tips to be removed. -1 for 2k");
    desc.AddOption("no_bubble", "", opt.no_bubble, "do not remove bubbles");
    desc.AddOption("min_standalone", "", opt.min_standalone, "min length of a standalone contig to output to final.contigs.fa");
    desc.AddOption("min_contig", "", opt.min_contig, "min length of contig to output");

    try {
        desc.Parse(argc, argv);

        if (opt.sdbg_name == "") {
            throw std::logic_error("no succinct de Bruijn graph name!");
        }
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: " << argv[0] << " -s sdbg_name -o output_prefix" << std::endl;
        std::cerr << "options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }
}

void PrintStat(Histgram<int64_t> &h) {
    // total length
    int64_t sum = h.sum();

    xlog("Total length: %lld, N50: %lld, Mean: %lld, number of contigs: %lld\n", (long long)sum, (long long)h.Nx(sum * 0.5), (long long)h.mean(), (long long)h.size());
    xlog("Maximum length: %llu\n", h.maximum());
}

int main_assemble(int argc, char **argv) {
    AutoMaxRssRecorder recorder;

    ParseAsmOption(argc, argv);

    SuccinctDBG dbg;
    xtimer_t timer;

    {
        // graph loading
        timer.reset();
        timer.start();
        xlog("Loading succinct de Bruijn graph: %s ", opt.sdbg_name.c_str());
        dbg.LoadFromMultiFile(opt.sdbg_name.c_str(), false);
        timer.stop();
        xlog_ext("Done. Time elapsed: %lf\n", timer.elapsed());
        xlog("Number of Edges: %lld; K value: %d\n", (long long)dbg.size, dbg.kmer_k);
    }

    {
        // set parameters
        if (opt.num_cpu_threads == 0) {
            opt.num_cpu_threads = omp_get_max_threads();
        }

        omp_set_num_threads(opt.num_cpu_threads);
        xlog("Number of CPU threads: %d\n", opt.num_cpu_threads);

        if (opt.max_tip_len == -1) {
            opt.max_tip_len = dbg.kmer_k * 2;
        }
    }

    if (opt.max_tip_len > 0) { // tips removal
        timer.reset();
        timer.start();
        assembly_algorithms::RemoveTips(dbg, opt.max_tip_len, opt.min_standalone);
        timer.stop();
        xlog("Tips removal done! Time elapsed(sec): %lf\n", timer.elapsed());
    }

    // remove bubbles
    if (!opt.no_bubble) {
        timer.reset();
        timer.start();
        uint32_t num_bubbles = assembly_algorithms::PopBubbles(dbg);
        timer.stop();
        xlog("Number of bubbles removed: %u/%u, Time elapsed(sec): %lf\n",
             num_bubbles, timer.elapsed());
    }

    // output contigs
    FILE *out_contig_file = OpenFileAndCheck(opt.contig_file().c_str(), "w");
    FILE *out_contig_info = OpenFileAndCheck((opt.contig_file() + ".info").c_str(), "w");
    // construct unitig graph
    timer.reset();
    timer.start();
    UnitigGraph unitig_graph(&dbg);
    Histgram<int64_t> hist;
    unitig_graph.InitFromSdBG(&hist, out_contig_file, opt.min_contig);
    timer.stop();
    xlog("unitig graph size: %u, time for building: %lf\n", unitig_graph.size(), timer.elapsed());


    PrintStat(hist);

    fprintf(out_contig_info, "%lld %lld\n", (long long)(hist.size()), (long long)(hist.sum()));

    fclose(out_contig_file);
    fclose(out_contig_info);

    return 0;
}