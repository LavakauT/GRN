# Informative GRN with kmer (Following KmerDiscovery, ML_pipeline, and kmer_similarity)




# Environment Requirements
- biopython 1.78
- matplotlib 3.5.3
- numpy 1.21.5
- pandas 1.3.5
- python 3.7.0
- scikit-learn 1.0.2
- scipy 1.7.3
```
$ wget http://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O ~/miniconda.sh
$ bash ~/miniconda.sh -b -p $HOME/miniconda
$ export PATH="$HOME/miniconda/bin:$PATH"
# source conda.sh
$ conda create -n ml python==3.7.0 
$ conda install biopython
$ conda install matplotlib
$ conda install pandas
$ conda install scikit-learn

```

# Basic ML Pipeline

Code provided to:

1. Clean your data (ML_preprocess.py)
2. Define a testing set (test_set.py)
3. Select the best subset of features to use as predictors (Feature_Selection.py)
4. Train and apply a classification (ML_classification.py) or regression (ML_regression.py) machine learning model
5. Assess the results of your model (output from the ML_classification/ML_regression scripts with additional options in scripts_PostAnalysis)


# Post ML
## Extract Machine learning results
Please follow [summary_ml_results.R](https://github.com/LavakauT/ML-pipeline/blob/main/summary_ml_results.R) to extract F1/AUC_ROC results in each groups and algorithms. Then, you can chose a best algorithm to extract the K-mers.

## PCC filter K-mers in 10 times training results
```
> suggested Folder arrangement
> Group
> -algorithm
> --imp files
```
First, please follow [get_kmer-imp_overlap.py](https://github.com/LavakauT/ML-pipeline/blob/main/get_kmer-imp_overlap.py)  to extract important K-mers among 10 times training in different algorithms.
Input: all imp files; Output: imp.txt (summary file including K-mers, imp value, and counts)


Then apply [ml_pcc_filte.R](https://github.com/LavakauT/ML-pipeline/blob/main/ml_pcc_filter.R) to pull out distinct and important K-mers in the selected algorithm.
Input: imp.txt; Output: imp_distinct_pcc_enriched_kmer.txt
