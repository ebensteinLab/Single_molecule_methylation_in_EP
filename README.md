# Single_molecule_methylation_in_EP:
Find methylation patterns of enhancer-promoter pairs coexisting on the same DNA molecule, from Bionano Genomics data.

### Preprocess:
Data: Optical methylation maps. https://genome.cshlp.org/content/29/4/646
#####Bionano outputs to single molecule methylation labels list:
Bionano .molecules file (.BNX) ->\
-> alignment (Bionano Access + Solve) to alignment files: .xmap, r.cmap, q.cmap.->\
-> alignment files to filtered single molecule methylation labels lists (.tsv): https://github.com/ebensteinLab/Irys-data-analysis \
-> single molecule methylation labels .tsv fields: \
molID, mol chr, mol start, mol end, a comma-separated list of labels locations in the query channel. 


###Single_molecule_methylation_in_EP:
The script 'Single_molecule_methylation_in_EP' looks for single DNA molecules that contain a full enhancer-promoter pair,
 and records the methylation states of the enhancer and promoter on the molecule 
 (elements that contain at least one label are marked "unmethylated", elements with not labels are marked "methylated").\
 Then the methylation states of each E-P pair from all the molecules containing it are summed. 
 

can process multiple single-molecule methylation labels lists from the same directory.
input files:
* directory of multiple single molecule methylation labels .tsv files
* promoters .BED file with the following fields (If processing ROM data, discard promoters that do no contain TCGA sequence sites)
chr prom_start  prom_end    gene_sym    other_field(TCGA count)
* enhancers file (If processing ROM data, discard enhancers that do no contain TCGA sequence sites):
chr enh_start  enh_end  gene_sym(linked to enhancer)    enh_ID    other_field(TCGA count)

other inputs:
* output directory
* minimal covergae (default = 10)
* minimal E-P distance (default = 5000bp)
* number of chromosomes to process (numerically sorted 1-24) (default =24)
```
usage: SM_methylation_combinations_series_of_inputs.py [-h] [-pr PROM]
                                                       [-en ENH]
                                                       [-m INPUT DIR]
                                                       [-o OUTPUT DIR]
                                                       [-n CHROMNUM]
                                                       [-mc MINCOV]
                                                       [-d DISTANCE]

get distributions of methylation combinations in enhancer-promoters pairs

optional arguments:
  -h, --help            show this help message and exit
  -pr PROM, --prom PROM
                        promoters file
  -en ENH, --enh ENH    enhancers file
  -m INPUT DIR, --molcs INPUT DIR
                        input files directory
  -o OUTPUT DIR, --output OUTPUT DIR
                        output files directory
  -n CHROMNUM, --chromnum CHROMNUM
                        number of chromosomes to process
  -mc MINCOV, --mincov MINCOV
                        minimal molecules number per pair
  -d DISTANCE, --distance DISTANCE
                        minimal distance between enh and prom [bp]
```


The final output for each input single-molecule labels list is a tab-separated txt file with the fields:
 Gene_sym, enh_ID, number of molecules with:\
  enh_m_prom_m, enh_um_prom_um, enh_um_prom_m, ehn_m_prom_um,
  chromosome, total number of molecules, and E-P distance.