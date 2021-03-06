## HM_and_TFbinding

These scripts can generate results in our manuscript entitled "Relationship between histone modifications and transcription factor binding is protein family specific". Here is the pipeline to model in vivo TF binding specificites:

1. Scan genome-wide transcription factors (TFs) binding sites (BSs) from ENCODE ChIP-seq data. Note that only the best hit per ChIP-sequence is considered;
2. Define non-BSs. Scan TF motifs from chromatin accessible regions, having exactly-matched core motif and distinct locations with BSs. Moreover, non-BSs have similar distribution of averaged chromatin accessibility with BSs.
3. Scan histone modification (HM) patterns around BSs and non-BSs in base pair resolution;
4. Build L2-regularized MLR models with different combinations of features.

## Dependencies

The pipeline requires:

* python 2.7 
* the BioPython module www.biopython.org
* BEDTools suite (Quinlan and Hall 2010): bedtools coverage (aka coverageBed, v2.17.0, a version that can take -abam option)
* the [FIMO package](http://meme-suite.org/doc/fimo.html) (Grant et al. 2011)
* the [DNAshapeR R package](http://bioconductor.org/packages/release/bioc/html/DNAshapeR.html) (Chiu et al. 2016) or the [DNAshape website](http://rohslab.cmb.usc.edu/DNAshape/) (Zhou et al. 2013)

## Tutorial

You can check MYC/MYC_HM_TFbinding_tutorial.ipynb to see how the pipeline is implemented for MYC. This pipeline allows us to obtain Figure 2 and 3, which are main results in our paper. Note that only TF binding data in the GM12878 cells are incldued here. If you are inteseted in generating other figures, or obtaining binding data in the K562 and H1-hESC cell lines, or reporting any problem, feel free to start an github issue or send an email to Beibei Xin at bxin@usc.edu.

## Project home page

For information on the source tree, examples, issues, and pull requests, see

    https://github.com/xinbeibei/HM_and_TFbinding/

## Cite

If you use any script or data from this GitHub folder, please cite:

Xin, B., & Rohs, R. (2018). Relationship between histone modifications and transcription factor binding is protein family specific. Genome research, 28(3), 321-333. doi: 10.1101/gr.220079.116. 
