#!/bin/bash

cd vars
python setup.py build_ext --inplace
cd ..
python setup.py build_ext --inplace

