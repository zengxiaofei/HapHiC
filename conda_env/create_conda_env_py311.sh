#!/bin/bash

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Creation Time: 2024-08-08 11:39



conda config --add channels conda-forge

# There is a serious memory leaking problem in Python 3.11.0
conda create -n haphic python=3.11.9
conda activate haphic

conda install mkl
conda install sparse_dot_mkl "numpy<2.0.0"

pip3 install scipy matplotlib
pip3 install scikit-learn networkx "pysam==0.20.0"

pip3 install portion

conda env export > environment_py311.yml
