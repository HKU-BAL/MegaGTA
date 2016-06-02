MegaGTA
=========

### Getting Started

```
git clone https://github.com/HKU-BAL/megagta.git
cd megagta
git submodule update --init --recursive
./make.sh

# assemble with default parameters
bin/megagta.py -r reads.fq -g gene_list.txt -o output_dir

# post-processing with clustering similarity threshold 0.99 (i.e. 1-0.01)
bin/post_proc.sh -g gene_resources_dir -d output_dir/contigs -m 16G -c 0.01

```

### Required Tools
* HMMER 3.1 (http://hmmer.janelia.org)
* UCHIME (http://drive5.com/usearch/manual/uchime_algo.html)
* ant 1.9.2 or greater (for compiling RDPTools)
* g++ 4.6 or greater (for compiling MegaGTA core)

### Simple Settings
* In `bin/post_proc.sh` users are supposed to modidy the paths for HMMER & UCHIME to their own ones.

* `gene_list.txt` is a config file. Each line of the file specifying interested genes and some of their information required by the assembly, including forward/reverse HMM, and the reference multiple sequence alignment of that gene family. Following is an example.
```
rplB share/RDPTools/Xander_assembler/gene_resource/rplB/rplB_forward.hmm share/RDPTools/Xander_assembler/gene_resource/rplB/rplB_reverse.hmm share/RDPTools/Xander_assembler/gene_resource/rplB/rplB_ref_alignment.faa
nirK share/RDPTools/Xander_assembler/gene_resource/nirK/nirK_forward.hmm share/RDPTools/Xander_assembler/gene_resource/nirK/nirK_reverse.hmm share/RDPTools/Xander_assembler/gene_resource/nirK/nirK_ref_alignment.faa
```

* `gene_resources_dir` is an input folder for post-processing which follows the format Xander assembler used. It can be found through:
```
megagta/share/RDPTools/Xander_assembler/gene_resource
```
* If genes not included are in interest, please follow the [README file of Xander][1] to prepare the resource.

[1]: https://github.com/rdpstaff/Xander_assembler/blob/master/README.md#per-gene-preparation-requires-biological-insight
