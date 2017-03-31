## HM_and_TFbinding

These scripts can generate results in our manuscript entitled "Relationship between histone modifications and transcription factor binding is protein family specific". Here is the pipeline to model in vivo TF binding specificites:

1. Scan genome-wide transcription factors (TFs) binding sites (BSs) from ENCODE ChIP-seq data. Note that only the best hit per ChIP-sequence is considered;
2. Define non-BSs. Scan TF motifs from chromatin accessible regions, having matched GC content and excluding genomic locations as BSs;
3. Scan histone modification (HM) patterns around BSs and non-BSs in base pair resolution;
4. Build L2-regularized MLR models with different combinations of features.

## Dependencies

The pipeline requires:

* python 2.7 
* BEDTools suite (Quinlan and Hall 2010)
* the [FIMO package](http://meme-suite.org/doc/fimo.html)
* the [BiasAway software](https://github.com/wassermanlab/BiasAway)
* the [DNAshapeR R package](http://bioconductor.org/packages/release/bioc/html/DNAshapeR.html) or the [DNAshape website](http://rohslab.cmb.usc.edu/DNAshape/)

## Tutorial

You can check MYC/MYC_HM_TFbinding_tutorial.ipynb to see how the pipeline is implemented for MYC. This pipeline allows us to obtain Figure 2 and 3, which are main results in our paper. If you are inteseted in generating other figures or you want to report any problems, feel free to start an github issue or send an email to Beibei Xin at bxin@usc.edu.

## Project home page

For information on the source tree, examples, issues, and pull requests, see

    https://github.com/xinbeibei/HM_and_TFbinding/
