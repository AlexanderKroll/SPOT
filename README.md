# Description
This repository contains the code and datasets to reproduce the results and figures and to train the models from our manuscript
 "A general substrate prediction model for transport proteins reveals differences between transporter types and localisations".


#### For people interested in using the trained prediction model, we implemented a [web server](https://spot.cs.hhu.de/) that allows an easy use of our trained model. The prediction tool can be run in a web-browser and does not require the installation of any software. Prediction results are usually ready within a few minutes.

#Requirements for running the code
I suggest to install a new conda environment with all required python packages:

#Creating new conda environment:
conda create -n Transporter python=3.7

#activate environment
conda activate Transporter

#Install jupyter notebook
conda install jupyter
jupyter notebook


#Installing python packages
conda create -n Transport python=3.7
conda activate Transport
conda install jupyter
conda install pandas
pip install goatools
pip install PubChemPy
pip install biopython
pip install openpyxl
pip install bioservices
conda install -c rdkit rdkit
pip install scikit-learn

conda install -c conda-forge py-xgboost=1.3.3
pip install torch


