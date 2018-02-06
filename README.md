# iFeature: A python package and web server for features extraction and selection from protein and peptide sequences

iFeature is a comprehensive Python-based toolkit for generating various numerical feature representation schemes from protein or peptide sequences. iFeature is capable of calculating and extracting a wide spectrum of 18 major sequence encoding schemes that encompass 53 different types of feature descriptors. Furthermore, iFeature also integrates five kinds of frequently used feature clustering algorithms, four feature selection algorithms and three dimensionality reduction algorithms. 
# Installation

  - Download iFeature by 
  ```sh
  git clone https://github.com/Superzchen/iFeature
  ```
  iFeature is an open-source Python-based toolkit, which operates depending on the Python environment (Python Version 3.0 or above) and can be run on multi-OS systems (such as Windows, Mac and Linux operating systems). Before running iFeature, user should make sure all the following packages are installed in their Python environment: sys, os, shutil, scipy, argparse, collections, platform, math, re, numpy (1.13.1), sklearn (0.19.1), matplotlib (2.1.0), and pandas (0.20.1). For convenience, we strongly recommended users to install the Anaconda Python 3.0 version (or above) in your local computer. The software can be freely downloaded from https://www.anaconda.com/download/.
# For users who want to generate descriptors by our provided iFeature package :
cd to the iFeature folder which contains iFeature.py, iFeaturePseKRAAC.py, cluster.py and feaSelector.py. All the functions regarding feature extraction, feature or sample clustering and feature selection analysis can be executed through these four main programs by specifying the parameter ‘--type’. 

“iFeature.py” is the main program used to extract 37 different types of feature descriptors:
```sh
python iFeature.py --help 
```
“iFeaturePseKRAAC.py” is the program used to extract the 16 types of pseudo K-tuple reduced amino acid composition (PseKRAAC) feature descriptors: 
```sh
python iFeaturePseKRAAC.py --help
``` 
Furthermore, the iFeature package contains other Python scripts to generate the position-specific scoring matrix (PSSM) profiles, predicted protein secondary structure and predicted protein disorder, which have also been often used to improve the prediction performance of machine learning-based classifiers in conjunction with sequence-derived information. The three dimensionality reduction algorithms are also included in the ‘scripts’ directory.
#### Examples for users to extract descriptors.
The input protein or peptide sequences for iFeature.py and iFeaturePseKRAAC.py should be in fasta format, Please find the input example in 'example' folder. Running the following command to obtain the `Amino Acid Composition (AAC)` descriptors:
```sh
python python iFeature.py --file examples/test-protein.txt --type AAC
```
'--file' specify the input file, while the '--type' is the descriptor type, the abbreviation of the descriptor types can be obtained by run `python iFeature.py --help` 

Users can generate different descriptors by change the descriptor type specified by '--type', For example, run the following command to generate the `Dipeptide Deviation from Expected Mean (DDE)` descriptor:
```sh
python python iFeature.py --file examples/test-protein.txt --type DDE
```
