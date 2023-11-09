#!/bin/bash

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Creation Time: 2022-09-14 21:58


conda create -n haphic python=3.10
conda activate haphic
pip3 install numpy scipy matplotlib
pip3 install sklearn networkx pysam

conda install -c intel mkl
conda config --add channels conda-forge
conda install sparse_dot_mkl

pip3 install portion

conda env export > environment_py310.yml
