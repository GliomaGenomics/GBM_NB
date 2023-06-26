# Analysis of nanobiopsy and whole-cell lysate sequencing datasets

## conda
```
conda activate env/nb3
```
## versions

STAR v2.7.3b
Other method versions are recorded in env/nb3.yml

## genome data
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz

TurboGFP added to end of GRCh38.primary_assembly.genome.fa

#scripts/collapse_annotation.py from https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py
python scripts/collapse_annotation.py downloaded_data/gencode/gencode.v27.annotation.gtf downloaded_data/gencode/gencode.v27.genes.gtf

## STAR
```
qsubsec scripts/star_index.qsubsec

# Whole-cell
for f in data/nano-biopsy_whole-cell/fastq/*L001_R1_001.fastq.gz ; do fi=$(basename $f _R1_001.fastq.gz) ; [ ! -f data/nano-biopsy_whole-cell/fastq/${fi}_R1_001_cutadapt.fastq.gz ] && qsubsec scripts/cutadapt.qsubsec DIR=data/nano-biopsy_whole-cell/fastq/ FILE=$fi -s ; done
for f in data/nano-biopsy_whole-cell/fastq/FACS_*_R1_*cutadapt.fastq.gz ; do fi=$(basename $f _R1_001_cutadapt.fastq.gz) ; qsubsec scripts/fastqc.qsubsec DIR=data/nano-biopsy_whole-cell/fastq/ FILE=$fi -s; done
qsubsec scripts/star_smartseq2.qsubsec SET=nano-biopsy_whole-cell NAME=1 MIN=0.66 
qsubsec scripts/bamsplit.qsubsec -s

# Nano-biopsy
for f in data/nano-biopsy/fastq/*L003_R1_001.fastq.gz ; do fi=$(basename $f _R1_001.fastq.gz) ; [ ! -f data/nano-biopsy/fastq/${fi}_R1_001_cutadapt.fastq.gz ] && qsubsec scripts/cutadapt.qsubsec DIR=data/nano-biopsy/fastq/ FILE=$fi -s ; done
qsubsec scripts/star_smartseq2.qsubsec NAME=1 MIN=0.66 SET=nano-biopsy -s
qsubsec scripts/bamsplit.qsubsec 
```

## Picard CollectRnaSeqMetrics
```
#setup
module load java/8u172 ; java -jar ~/picard.jar CreateSequenceDictionary R=downloaded_data/gencode/GRCh38.primary_assembly.genome.fa O=downloaded_data/gencode/GRCh38.primary_assembly.genome.dict
grep --color=none -i rrna downloaded_data/gencode/gencode.v27.annotation.gtf > downloaded_data/gencode/gencode.v27.annotation_rib.gtf
export PATH="./bin:$PATH" ; gff2bed < downloaded_data/gencode/gencode.v27.annotation_rib.gtf > downloaded_data/gencode/gencode.v27.annotation_rib.bed
module load java/8u172 ; java -jar ~/picard.jar BedToIntervalList I=downloaded_data/gencode/gencode.v27.annotation_rib.bed O=downloaded_data/gencode/gencode.v27.annotation_rib.interval_list SD=downloaded_data/gencode/GRCh38.primary_assembly.genome.dict
./gtfToGenePred -genePredExt downloaded_data/gencode/gencode.v27.annotation.gtf downloaded_data/gencode/gencode.v27.annotation.ref_flat.txt
#cat downloaded_data/gencode/gencode.v27.annotation.ref_flat.txt | awk '{print $12"\t"$0}' | cut -d$'\t' -f1-11 > downloaded_data/gencode/temp.txt ; mv downloaded_data/gencode/temp.txt downloaded_data/gencode/gencode.v27.annotation.ref_flat.tximport

#run
DIR=nano-biopsy/star_gfp_min0.66_minNmatch3
for f in data/${DIR}/split_bams/Aligned.out_B10_180122_M059K_S202_L00* ; do fi=$(basename $f) ; qsubsec scripts/picard_rnaseqstats.qsubsec DIR=data/${DIR}/split_bams FILE=$fi  ; done
for f in data/${DIR}/split_bams/*stats.txt ; do head -n 8 $f | tail -n 2 > ${f}_g ;done
python scripts/merge.py -i data/${DIR}/split_bams/ -o data/${DIR}/bam_stats/ -s _stats.txt_g

for f in data/${DIR}/split_bams/*unique.bam ; do echo $(basename $f _sorted_rg_unique.bam) ; done
for f in data/${DIR}/split_bams/*rg_flagstat* ; do grep 'primary mapped' $f | cut -f1 -d '+' ; done
for f in data/${DIR}/split_bams/*rg_unique_flagstat* ; do grep 'paired in' $f | cut -f1 -d '+' ; done
for f in data/${DIR}/split_bams/*rg_flagstat* ; do grep 'paired in' $f | cut -f1 -d '+' ; done
for f in data/${DIR}/split_bams/*unique.bam ; do samtools view $f 'TurboGFP' | wc -l ; done
```

## seurat
```
scripts/seurat.R
#protein coding genes taken from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.pc_translations.fa.gz
grep '>' downloaded_data/gencode/gencode.v27.pc_translations.fa | tr '|' '\t' | cut -f 3 | sort -u > downloaded_data/gencode/protein_coding_v27.txt
```

## GSVA

Gene sets for mesenchymal and proneural phenotypes were taken from the top and bottom 50 weighted genes from PC1 for the 10x single-cell dataset in Wang et al. 2019 (https://pubmed.ncbi.nlm.nih.gov/31554641/).
The stemness gene set was taken from Patel et al. 2014 (https://pubmed.ncbi.nlm.nih.gov/24925914/).
```
#run
gsva.R

### manually process results into paired format in 'wang/gsva_paired_results.txt'

#plot paired
gsva_plot.R
```

##Figures
Figures for the manuscript are in the wang and seurat folders. 
The QC filtering diagram is in meta/.
