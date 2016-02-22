megaGTA
=========

## Getting Started

```
git clone git@bitbucket.org:almighty_mimicry/kingassembler.git
cd kingAssembler
git submodule update --init --recursive
./make.sh

# assemble with default parameters
bin/megagta.py -r reads.fq -g gene_list.txt -o output_dir

# post-processing with clustering similarity threshold 0.99 (i.e. 1-0.01)
bin/post_proc.sh -g gene_resources_dir -d output_dir/contigs -m 16G -c 0.01

```

* `gene_list.txt` is a config file for specifying interested genes and some of their information required by the assembly:
		rplB rplB_forward.hmm rplB_reverse.hmm rplB_ref_alignment.faa
		nirK nirK_forward.hmm nirK_reverse.hmm nirK_ref_alignment.faa
		...

* `gene_resources_dir` is an input folder for post-processing which follows the format Xander assembler used. It can be found through:
		megahit_gt/share/RDPTools/Xander_assembler/gene_resource
If genes not included are in interest, please follow the [README][1] file of Xander to prepare the resource.

[1]: https://github.com/rdpstaff/Xander_assembler/blob/master/README.md
