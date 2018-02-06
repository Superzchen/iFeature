# iFeature: A python package and web server for features extraction and selection from protein and peptide sequences

iFeature is a comprehensive Python-based toolkit for generating various numerical feature representation schemes from protein or peptide sequences. iFeature is capable of calculating and extracting a wide spectrum of 18 major sequence encoding schemes that encompass 53 different types of feature descriptors. Furthermore, iFeature also integrates five kinds of frequently used feature clustering algorithms, four feature selection algorithms and three dimensionality reduction algorithms. 
# Installation

  - Download iFeature by 
  ```sh
  git clone https://github.com/Superzchen/iFeature
  ```
  iFeature is an open-source Python-based toolkit, which operates depending on the Python environment (Python Version 3.0 or above) and can be run on multi-OS systems (such as Windows, Mac and Linux operating systems). Before running iFeature, user should make sure all the following packages are installed in their Python environment: sys, os, shutil, scipy, argparse, collections, platform, math, re, numpy (1.13.1), sklearn (0.19.1), matplotlib (2.1.0), and pandas (0.20.1). For convenience, we strongly recommended users to install the Anaconda Python 3.0 version (or above) in your local computer. The software can be freely downloaded from https://www.anaconda.com/download/.
# For general users who want to generate descriptors by our provided iFeature package :
cd to the iFeature folder which contains iFeature.py, iFeaturePseKRAAC.py, cluster.py and feaSelector.py. All the functions regarding feature extraction, feature or sample clustering and feature selection analysis can be executed through these four main programs by specifying the parameter '--type'. 

"iFeature.py" is the main program used to extract 37 different types of feature descriptors. For details of other parameters, run:
```sh
python iFeature.py --help 
```
"iFeaturePseKRAAC.py" is the program used to extract the 16 types of pseudo K-tuple reduced amino acid composition (PseKRAAC) feature descriptors. For details of other parameters, run: 
```sh
python iFeaturePseKRAAC.py --help
```
"cluster.py" is the program used for running the feature or sample clustering algorithms. For details of other parameters, run:
```sh
python cluster.py --help
```
"feaSelector.py" is the fourth main program used to implement the feature selection algorithms. For details of other parameters, run:
```sh
python feaSelector.py --help
```
Furthermore, the iFeature package contains other Python scripts to generate the position-specific scoring matrix (PSSM) profiles, predicted protein secondary structure and predicted protein disorder, which have also been often used to improve the prediction performance of machine learning-based classifiers in conjunction with sequence-derived information. The three dimensionality reduction algorithms are also included in the `scripts` directory.
#### Examples for users to extract descriptors from `iFeature.py`. All files in the example commands can be found in the `examples` directory. 
The input protein or peptide sequences for iFeature.py and iFeaturePseKRAAC.py should be in fasta format, Please find the example in `example` folder. The following parameters are required by `iFeature.py`:
* --help    show help of 'iFeature.py'
* --file    protein/peptide sequence file in fasta format
* --type    feature types for protein sequence analysis
* --path    data file path used for 'PSSM', 'SSEB(C)', 'Disorder(BC)', 'ASA' and 'TA' encodings
* --train   training file in fasta format only used for 'KNNprotein' or 'KNNpeptide' encodings
* --label   sample label file only used for 'KNNprotein' or 'KNNpeptide' encodings
* --order   output order for of Amino Acid Composition (i.e. AAC, EAAC, CKSAAP, DPC, DDE, TPC) descriptors
* --userDefinedOrder user defined output order for of Amino Acid Composition (i.e. AAC, EAAC, CKSAAP, DPC, DDE, TPC) descriptors
* --out     the generated descriptor file
Running the following command to obtain the `Composition of k-spaced Amino Acid Pairs (CKSAAP)` descriptor:
```sh
python python iFeature.py --file examples/test-protein.txt --type CKSAAP
```
Generally, users can generate different descriptors by change the descriptor type specified by '--type', For example, run the following command to generate the `Dipeptide Deviation from Expected Mean (DDE)` descriptor:
```sh
python python iFeature.py --file examples/test-protein.txt --type DDE
```
For some descriptors (e.g. `PSSM`, `Disorder`, `ASA`, `TA` and `SSEB`), the predicted protein property file should be supplied by the `--path` parameter. For example, run the following command to generate `PSSM` descriptor:
```sh
python iFeature.py --file examples/test-peptide.txt --type PSSM --path examples/predictedProteinProperty
``` 
`KNNprotein` and `KNNpeptide` descriptors requires an extra training file and a label file, which is spedified by `--train` and `--label`. Run the following command to generate the `KNNprotein` descriptor:
```sh
python iFeature.py --file examples/test-peptide.txt --type KNNpeptide --train examples/train-peptide.txt --label examples/label.txt
``` 
For the six descriptors in `Amino Acid Composition` group, user can specify the output order by `--order` and `--userDefinedOrder`, three amino acids order (i.e. alphabetically, polarity and side chaim volume) were supplied by iFeature. Run the following command to generate the AAC descriptor with the 'polarity' order:
```sh
python iFeature.py --file examples/test-protein.txt --type AAC --order polarity
```
Run the following command to generate the `AAC` descriptor with a user-defined order:
```sh
python iFeature.py --file examples/test-protein.txt --type AAC --order userDefined --userDefinedOrder YWVTSRQPNMLKIHGFEDCA
```
For some of the descriptors, user can adjust the default parameters by advanced usage. The detatiled advanced usage for each type of descriptor can be found in `iFeature_manual.pdf` in `iFeature` directory. For example, for `CKSAAP` descriptor, advanced users can adjust the size of the sliding window to <N> (the default is 5) by running the following Python command:
```sh
python codes/EAAC.py examples/test-peptide.txt 3 EAAC.tsv
```
The default output file is `encoding.tsv`, which can be specified by `--out`. For example:
```sh
python iFeature.py --file examples/test-protein.txt --type AAC --out AAC.txt
```
#### Examples for users to extract descriptors from `iFeaturePseKRAAC.py`.
The 16 types of reduced amino acid alphabets with different clustering approaches can be used to generate different versions of pseudo reduced amino acid compositions (PseRAACs).

The following parameters are required by `iFeaturePseKRAAC.py`:
* --file protein/peptide sequence file in fasta format
* --type descriptor type in the PseKRAAC feature group
* feature types for protein sequence analysis, two alternative modes (g-gap and lambda-correlation) are available, with the ‘g-gap’ model as the default. 
* K-tuple value, three K-tuple values (i.e. 1, 2 and 3) are available, default is 2
* --gap_lambda gap value for the `g-gap` model or lambda value for the `lambda-correlation’ model, 10 values are available (i.e. 0, 1, 2, …, 9)
* the reduced amino acids cluster type.
Users can run the following command to view the available values for each descriptor type:
```sh
python iFeaturePseKRAAC.py --show
```
Use the following command to extract the PseKRAAC feature descriptors:
```sh
python iFeaturePseKRAAC.py --file examples/test-protein.txt --type type1 --subtype lambda-correlation --ktuple 2 --gap_lambda 2 --raactype 5
```





