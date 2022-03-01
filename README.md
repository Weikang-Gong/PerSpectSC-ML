# **Persistent spectral simplicial complex-based machine learning (PerSpectSC-ML) reveals principles of chromosomal structure in cellular differentiation**
This manual is for the code implementation of the paper "Persistent spectral simplicial complex-based machine learning (PerSpectSC-ML) reveals principles of chromosomal structure in cellular differentiation".

## Preparation
Code Requirements
Platform: Python>=3.6, MATLAB 2019B
Python Packages needed: math, numpy>=1.19.5, scipy>=1.4.1, scikit-learn>=0.20.3, GUDHI 3.0.0

## Details about each step

**Step 1:** Distance matrix
Before the representation, the distance matrix based on Hi-C data is got through HiC_TDA.py (https://github.com/Kingsford-Group/hictda).

**Step 2:** Persistent Spectral simplicial complex(PerSpectSC) representation and Feature generation
For each chromosome, the distance matrix is used to construct the simplicial complexes to generate the Hodge Laplacians.

Chromosome_Hodge_Laplacian_L0.py is used to compute the eigenvalues of a 0-dimensional Laplacian matrix.
Chromosome_Hodge_Laplacian_L1.py is used to compute the eigenvalues of a 1-dimensional Laplacian matrix.

Persistent_Attributes_Structural_Classification.m is used for chromosomal descriptors generation.

**Step 3:** t-SNE-assisted k-means for classification of different cell types
Persistent_Attributes_Structural_Classification.m is used for the classification of different cell types.

## Help
For any questions, please contact us by weikanggong@emails.bjut.edu.cn.