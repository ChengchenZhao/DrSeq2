# DrSeq2.0

Quality control and analysis pipeline for parallel single cell transcriptome and epigenome data

By applying this pipeline, DrSeq2 take sequencing files as input and provides four groups of QC measurements for given data, including reads level, bulk-cell level, individual-cell level and cell-clustering level QC.
Here we provide an example to get you easily started on a linux/MacOS system with python and R installed. To run DrSeq2 with options specific to your data, you need to see Manual section for detailed usage.

Step1.Install pipeline
-------------------------------------------------------------------------------------------------------------------------------------------
1. Make sure you have python2.7 and R(version >= 2.14.1) on linux or MAC OSX environment.
2. Get Dr.seq 
```shell
git clone https://github.com/ChengchenZhao/DrSeq2.0.git
```
3. Install Dr.seq on your server/computer
```shell
cd DrSeq2.0/bx-python-0.7.1
python setup.py install --prefix ~/bin
cd ../
python setup.py install --prefix ~/bin
```
