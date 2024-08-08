#!/bin/bash

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Creation Time: 2024-08-08 11:54

conda config --add channels conda-forge

conda create -n haphic_py312 python=3.12.4
conda activate haphic_py312

conda install mkl

conda install sparse_dot_mkl "numpy<2.0.0"

# pip3 install "numpy<2.0"
pip3 install scipy matplotlib

pip3 install scikit-learn networkx

# higher version of GCC may be required
pip3 install pysam

pip3 install portion

conda env export > environment_py312.yml
