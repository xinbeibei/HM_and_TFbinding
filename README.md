## HM_and_TFbinding

This pipeline allows for modeling in vivo TF binding specificites:

1. Scan genome-wide transcription factors (TFs) binding sites (BSs) from ENCODE ChIP-seq data;
2. Define non-BSs. Scan TF motifs from chromatin accessible regions, having matched GC content and excluding genomic locations as BSs;
3. Scan histone modification (HM) patterns around BSs and non-BSs in base pair resolution;
4. Build L2-regularized MLR models with different combinations of features.

Note that only the best hit per ChIP-sequence is considered in the current version of the module.

## Dependencies

The pipeline requires:

* python2.7 
* the BioPython module www.biopython.org.
* the [TFFM package](http://cisreg.cmmt.ubc.ca/TFFM/doc/index.html) accessed from your
PYTHONPATH environment variable.
* the [scikit-learn module](http://scikit-learn.org/stable).
* the [pandas module](http://pandas.pydata.org).
* access to bigWig files providing the values of the DNA shape features HelT,
MGW, ProT, and Roll from your genome interest along with the second order
computation these features. Please visit the
[GBshape website](rohsdb.cmb.usc.edu/GBshape).
* the [bwtool](https://github.com/CRG-Barcelona/bwtool).

## Tutorial

You can find some examples of how to run the DNAshapedTFBS.py tool in the script
test.sh provided in the test/ repository of this package.

The script feature_importance_heatmap.py plots the heatmap(s) of trained
classifier(s). Note that the current version only works for PSSM/TFFM + DNA
shape classifiers. You can get help on how to use it by typing

python2.7 feature_importance_heatmap.py -h

## Project home page

For information on the source tree, examples, issues, and pull requests, see

    http://github.com/amathelier/DNAshapedTFBS

## Cite

If you use the DNAshapedTFBS tool, please cite

* A. Mathelier, B. Xin, T.-P. Chiu, L. Yang, R.R. Rohs, and W.W. Wasserman (2016)
DNA shape features improve transcription factor binding site predictions *in
vivo*. *Cell Systems*,
DOI:[10.1016/j.cels.2016.07.001](http://dx.doi.org/10.1016/j.cels.2016.07.001).
