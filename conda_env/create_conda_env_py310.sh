#!/bin/bash

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Creation Time: 2024-08-08 11:10


conda config --add channels conda-forge

conda create -n haphic python=3.10.14
conda activate haphic

conda install mkl

conda install sparse_dot_mkl "numpy<2.0.0"

pip3 install scipy matplotlib
pip3 install scikit-learn networkx "pysam==0.20.0"

pip3 install portion

conda env export > environment_py310.yml
