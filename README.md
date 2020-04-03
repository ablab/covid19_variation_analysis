This repo contains support scripts for biorxiv draft LINK.

Here you can find three python scripts:

construct_vcf_from_multifasta.py - construct vcf file from fasta file with complete SARS-CoV-2 sequences downloaded from [GISAID](https://www.gisaid.org/) database  

check_singletons_db.py - script to compare whether variations in your sample were already seen in GISAID db (requires vcf constructed by previous script)  

find_singleton_mutations.py - finds all singleton mutations in vcf file constructed by first script and output only those ones that were sequenced before given date.  


Usage and parameters are available by running these scripts without any options.  


First two scripts requires [minimap2](https://github.com/lh3/minimap2/) and [picard](https://broadinstitute.github.io/picard/) to be either in $PATH or provided as an option.
Also, these scripts requires k8 to run paftools.js (which comes within minimap2 distribution, see details [here](https://github.com/lh3/minimap2/blob/master/misc/README.md) )

NB! picard merges vcf containing multiple variations in the same place and different samples into format   
MN908947.3      13     .       T       G       60      .       SAMPLE1       GT      1/1  
MN908947.3      13     .       T       G       60      .       SAMPLE2       GT      1/1  
MN908947.3      13     .       T       G       60      .       SAMPLE3       GT      1/1  
This is not standard for vcf format, so we do not guarantee that such vcfs would be suitable for any external vcf-processing tools.


Also, there is available a draw_figure_1.ipynb jupyter notebook that can be used to rebuilt figure 1 from the draft.
