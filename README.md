## *iLearn* - the updated version of *iFeature* is now available! 
*iLearn* is a Python Toolkit and Web Server Integrating the Functionality of Feature Calculation, Extraction, Clustering, Feature Selection, Feature Normalization, Dimension Reduction and Model Construction for Classification, Best Model Selection, Ensemble Learning and Result Visualization for DNA, RNA and Protein Sequences. Please refer to https://github.com/Superzchen/iLearn for details.

# *iFeature*: A python package and web server for features extraction and selection from protein and peptide sequences

*iFeature* is a comprehensive Python-based toolkit for generating various numerical feature representation schemes from protein or peptide sequences. *iFeature* is capable of calculating and extracting a wide spectrum of 18 major sequence encoding schemes that encompass 53 different types of feature descriptors. Furthermore, *iFeature* also integrates five kinds of frequently used feature clustering algorithms, four feature selection algorithms and three dimensionality reduction algorithms. 
# Installation

  - Download *iFeature* by 
  ```sh
  git clone https://github.com/Superzchen/iFeature
  ```
  *iFeature* is an open-source Python-based toolkit, which operates depending on the Python environment (Python Version 3.0 or above) and can be run on multi-OS systems (such as Windows, Mac and Linux operating systems). Before running *iFeature*, user should make sure all the following packages are installed in their Python environment: sys, os, shutil, scipy, argparse, collections, platform, math, re, numpy (1.13.1), sklearn (0.19.1), matplotlib (2.1.0), and pandas (0.20.1). For convenience, we strongly recommended users to install the Anaconda Python 3.0 version (or above) in your local computer. The software can be freely downloaded from https://www.anaconda.com/download/.
# For users who want to generate descriptors by our provided *iFeature* package :
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
Furthermore, the *iFeature* package contains other Python scripts to generate the position-specific scoring matrix (PSSM) profiles, predicted protein secondary structure and predicted protein disorder, which have also been often used to improve the prediction performance of machine learning-based classifiers in conjunction with sequence-derived information. The three dimensionality reduction algorithms are also included in the `scripts` directory.
### Examples for users to extract descriptors from `iFeature.py`. All files in the example commands can be found in the `examples` directory. 
The input protein or peptide sequences for iFeature.py and iFeaturePseKRAAC.py should be in fasta format, Please find the example in `example` folder. The following parameters are required by `iFeature.py`:
* `--help`    show help of 'iFeature.py'
* `--file`   protein/peptide sequence file in fasta format
* `--type`    feature types for protein sequence analysis
* `--path`    data file path used for 'PSSM', 'SSEB(C)', 'Disorder(BC)', 'ASA' and 'TA' encodings
* `--train`   training file in fasta format only used for 'KNNprotein' or 'KNNpeptide' encodings
* `--label`   sample label file only used for 'KNNprotein' or 'KNNpeptide' encodings
* `--order`   output order for of Amino Acid Composition (i.e. AAC, EAAC, CKSAAP, DPC, DDE, TPC) descriptors
* `--userDefinedOrder` user defined output order for of Amino Acid Composition (i.e. AAC, EAAC, CKSAAP, DPC, DDE, TPC) descriptors
* `--out`     the generated descriptor file
Running the following command to obtain the `Composition of k-spaced Amino Acid Pairs (CKSAAP)` descriptor:
```sh
python python iFeature.py --file examples/test-protein.txt --type CKSAAP
```
Generally, users can generate different descriptors by changing the descriptor type by specifying '--type', For example, run the following command to generate the `Dipeptide Deviation from Expected Mean (DDE)` descriptor:
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
For the six descriptors in `Amino Acid Composition` group, user can specify the output order by `--order` and `--userDefinedOrder`, three amino acids order (i.e. alphabetically, polarity and side chaim volume) were supplied by *iFeature*. Run the following command to generate the AAC descriptor with the 'polarity' order:
```sh
python iFeature.py --file examples/test-protein.txt --type AAC --order polarity
```
Run the following command to generate the `AAC` descriptor with a user-defined order:
```sh
python iFeature.py --file examples/test-protein.txt --type AAC --order userDefined --userDefinedOrder YWVTSRQPNMLKIHGFEDCA
```
For some of the descriptors, user can adjust the default parameters by advanced usage. The detatiled advanced usage for each type of descriptor can be found in `iFeatureManual.pdf` in `iFeature` directory. For example, for `CKSAAP` descriptor, advanced users can adjust the size of the sliding window to 3 (the default is 5) by running the following Python command:
```sh
python codes/EAAC.py examples/test-peptide.txt 3 EAAC.tsv
```
The default output file is `encoding.tsv`, which can be specified by `--out`. For example:
```sh
python iFeature.py --file examples/test-protein.txt --type AAC --out AAC.txt
```
### Examples for users to extract descriptors from `iFeaturePseKRAAC.py`.
The 16 types of reduced amino acid alphabets with different clustering approaches can be used to generate different versions of pseudo reduced amino acid compositions (PseRAACs).

The following parameters are required by `iFeaturePseKRAAC.py`:
* `--file` protein/peptide sequence file in fasta format
* `--type` descriptor type in the PseKRAAC feature group
* `feature` types for protein sequence analysis, two alternative modes (g-gap and lambda-correlation) are available, with the ‘g-gap’ model as the default. 
* `K-tuple` value, three K-tuple values (i.e. 1, 2 and 3) are available, default is 2
* `--gap_lambda` gap value for the `g-gap` model or lambda value for the `lambda-correlation’ model, 10 values are available (i.e. 0, 1, 2, …, 9)
* the reduced amino acids cluster type.
Users can run the following command to view the available values for each descriptor type:
```sh
python iFeaturePseKRAAC.py --show
```
Use the following command to extract the PseKRAAC feature descriptors:
```sh
python iFeaturePseKRAAC.py --file examples/test-protein.txt --type type1 --subtype lambda-correlation --ktuple 2 --gap_lambda 2 --raactype 5
```
### Examples for Feature analysis.
*iFeature* integrates several commonly used and useful clustering, feature selection and dimensionality reduction algorithms. In order to facilitate the understanding of the results for non-experts, a scatter diagram will be ploted to show the distribution of the cluster result. When ploting the scatter diagram, the t-SNE algorithm was used to reduce the high-dimensional descriptors to 2 dimensions, and then the scatter diagram was plotted by using the reduced descriptors and clustering result. The clustering result will be also stored in a file with text format. The file contained the information of clusters number and ownership for each sample. The clustering algorithms can be run by the following command:
```sh
python cluster.py --file descriptor.tsv --type <clustering_algorithm> --sof <sample/feature>
```
* `--file` is the descriptor output file, which is generated by `iFeature.py` or `iFeaturePseKRAAC.py`
* `--type` specifies the clustering algorithm (kmeans, hcluster, apc, meanshift and dbscan)
* `--sof` specifies the option for performing clustering for samples or features (default: samples)
Use the following command to perform the K-Means clustering:
```sh
python cluster.py --file examples/example.tsv --type kmeans --sof feature --nclusters 2
```
The feature selection algorithms can be run by the following command:
```sh
python feaSelector.py --file descriptor.tsv --type <feature_selection_algorithm> --label sample_label_file
```
* `--label` specify the sample class of the samples in `descriptor.tsv`.
Use the following command to perform the Chi-Square feature selection:
```sh
python feaSelector.py --file examples/example.tsv --type CHI2 --label examples/label.txt --out CHI2_feature.txt
```
After running the feature selection algorithm  the descriptors will be ranked according to their importance, the higher the ranking is, the more important the descriptor is.

In addition, three dimensionality reduction algorithms (PCA,LDA and t-SNE) have been implemented in the *iFeature* package and can be run using the following commands. To facilitate understanding, the dimensionality reduction result can be
visualized by the scatter diagram.

Use the following command to perform the PCA analysis:
```sh
python scripts/pcaAnalysis.py --file examples/example.tsv --ncomponents 3 --out pcaResult.txt
```
* `--ncomponents` specifies the number of principal components, `--out` specifies the name of the output file of PCA.
Use the following command to perform the LDA analysis:
```sh
python scripts/ldaAnalysis.py --file examples/example.tsv --ncomponents 3 --label examples/label.txt --out ldaResult.txt
```
Use the following command to perform the t-SNE analysis:
```sh
python scripts/tsneAnalysis.py --file examples/example.tsv
```

### Citation：
If you find *iFeature* useful, please kindly cite the following paper:

Chen Z, Zhao P, Li F, Leier A, Marquez-Lago TT, Wang Y, Webb GI, Smith AI, Daly RJ* , Chou KC* , Song J* . *iFeature*: a python package and web server for features extraction and selection from protein and peptide sequences. Bioinformatics, 2018, doi: 10.1093/bioinformatics/bty140.
