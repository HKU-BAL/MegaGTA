MEGAHIT-GT
=========

## Getting Started

```
git clone git@bitbucket.org:almighty_mimicry/kingassembler.git
cd kingAssembler
git submodule update --init --recursive
./make.sh
bin/prepare_gene_ref.sh gene_resource.txt
python bin/kingassembler.py --k-list 29,44 -r reads.fa --min-count 1 --gene_list gene_resource/gene_list.txt -o result
bin/megahit_gt_post_proc.sh -d result -h 16G -c 0.01

```

