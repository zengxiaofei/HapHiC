#!/bin/bash

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Creation Time: 2023-04-04 19:53


# There is a serious memory leaking problem in Python 3.11.0
conda create -n haphic_py311 python=3.11.2
conda activate haphic_py311

pip3 install numpy scipy matplotlib
pip3 install scikit-learn networkx pysam

conda install -c intel mkl
conda config --add channels conda-forge
conda install sparse_dot_mkl

pip3 install portion

conda env export > environment_py311.yml
