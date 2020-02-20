paper review
================
Yuan Xia
19/02/2020

### Goals and findings of the paper.

Researchers suspect that hereditary deafness is often caused by
pathogenic variants in hair cell genes. Subsequently, they have
conducted a RNA-Seq study which explored the gene expression in hair
cell during mouse inner ear development. different gene expression
patterns are found in different cell types (hair cells and surrounding
cells), tissues sourses (cochlea and utricle) and development
stages(E16, P0, P4 and P7). Additionally, the expression of genes
critical for mechanotransduction and deafness are founded in hair cell
which shows the possibility of the discovery of unknown deafness genes
in hair cells.

### Datasets, experiments and relevant analysis

1.  Since the hair cells are embedded in sensory epithelium with a
    variety of other cell types, enzymatic treatment was used to
    dissociate hair cells from surrounding cells. The FACs technique was
    used to further purify hair cells using the brighest GFP
    fluorescence signal to collect hair cells and the lowest
    fluorescence signal to collect other cells.

2.  The dataset contains the expression of 20,207 REfSeq mouse genes
    from 16 samples with three factors: cell types, tissue sources and
    developmental stages. The sequencing libraries converted from
    3’-enriched transcriptome tags that have similar gene length. Gene
    expression levels were measured using RPKM (Reads Per kilobase of
    thanscript per Milion mapped reads) in log2 scale.

3.  Several data visulization methods were used to find interesting
    differences between factor levels. FIrstly, Principle component
    analysis was used to plot the normalized gene expression data though
    reduceing the number of axes needed to display the important aspects
    of the data. Secondly, DESeq package was also used to identify
    differentially expressed genes between hair and surrounding cells
    and cochlear and utricular hair cells. Finally, histograms were used
    to display and compare the normalized number of reads for sepcific
    genes in different cells, tissue sources and devlopment stages.

4.  Fold change(FC) value was used to represent how big the relative
    difference is between selected highly enriched hair (GFP+) and
    surrounding cells (GFP-). It is defined as the GFP+/GFP- counts
    ratio.

5.  Flase discovery rate(FDR) was determined using the Benjamini and
    Hochberg method. It shows the rate of type I errors in hypothesis
    testing when conducting mutiple comparsions between different gene
    expression across cells, organs and ages.

6.  Hierarchical clustering heat maps were used to show patial and
    temporal expression patterns of genes at various developmetal
    stages. The distance of different genes are defined as correlation
    of two genes across samples, and the centroid linkage method was
    used to determine the distance between gene clusters.

### Comments and discussions

This study intented to comprehend the hair cells gene expression
specific for vestibular and cochlear HCs during mouse development and to
discover unknown deafness genes. Serval pairwise and multiple
comparsions were conducted by ploting the data to illustrate different
gene expression patterns. However, I cannot find their differential
expression analysis models or any p-values to show statistical
significance. The analytical mothods are rather simple, and most of
analysis steps are still at the stage of exploratory data analysis.
Statistical models are suggested to be used to further confirm the
findings.