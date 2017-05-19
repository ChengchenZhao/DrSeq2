# DrSeq2.0

Quality control and analysis pipeline for parallel single cell transcriptome and epigenome data

By applying this pipeline, DrSeq2 takes sequencing files as input and provides four groups of QC measurements for given data, including reads level, bulk-cell level, individual-cell level and cell-clustering level QC.
Here we provide an example to get you easily started on a linux/MacOS system with python and R installed. To run DrSeq2 with options specific to your data, you need to see Manual section for detailed usage.

STEP1.Install pipeline
-------------------------------------------------------------------------------------------------------------------------------------------
1. Make sure you have Python (version 2.7) and R (version 2.14.1 or higher) in a Linux or MacOS environment.
2. Get Dr.seq 
```shell
git clone https://github.com/ChengchenZhao/DrSeq2.0.git
```
3. Install Dr.seq on your server/computer
```shell
cd DrSeq2.0
```
for root user
```shell
sudo python setup.py install
```
If you are not a root user, DrSeq2 can be installed at a specific location with write permission.
```shell
$ python setup.py install --prefix /home/DrSeq2    # here you can replace "/home/DrSeq2" with any location you want
$ export PATH=/home/DrSeq2/bin:$PATH    # setup PATH, so that system knows where to find executable files
$ export PYTHONPATH=/home/DrSeq2/lib/python2.7/site-packages:$PYTHONPATH    # setup PYTHONPATH, so that DrSeq2 knows where to import modules
```
NOTE: To install DrSeq2 on MacOSX, users need to download and install Xcode beforehand.
DrSeq2 may take several minutes for installation.
Type:
```shell
$ DrSeq --help
```
If you see help manual, you have successfully installed DrSeq2.

STEP2.Prepare annotation and required softwares
-------------------------------------------------------------------------------------------------------------------------------------------
1. Obtain gene annotation information for the species of your DrChIP data or ATAC-seq data. 
We have provided a convenient link of pre-compiled gene annotation tables with genome background annotations for human (hg18) and mouse (mm9) data on our webpage.You can download annotation tables for other species from CEAS.(Skip 2 and 3 if you already have bowtie2 and bowtie2 index)
(Skip 4 if you already have macs14 installed)
2. Prepare the mapping software (We use bowtie2 for the quick start mode). 
We have provided a convenient link to an executable version of bowtie2 for Linux and MacOS users at out webpage (otherwise you must download the full bowtie2 package and compile it yourself, see the Manual section): 

For root users, simply copy executable bowtie2 to any default PATH (for example: /usr/local/bin).
```shell
$ unzip bowtie2-2.2.6-linux-x86_64.zip    # bowtie2-2.2.6-macos-x86_64.zip for macOS
$ cd bowtie2-2.2.6-linux-x86_64    # bowtie2-2.2.6-macos-x86_64 for macOS
$ sudo cp bowtie2* /usr/local/bin    # don’t forget to type * mark here
```
If you are not a root user, you can copy bowtie2 to the PATH (/home/DrSeq2/bin) you setup in step 1.
```shell
$ unzip bowtie2—2.2.6-linux-x86_64.zip    # bowtie2-2.2.6-macos-x86_64.zip for macOS
$ cd bowtie2-2.2.6-linux-x86_64    # bowtie2-2.2.6-macos-x86_64 for macOS
$ cp bowtie2* /home/DrSeq2/bin    # don’t forget to type * mark here
```
Type

```shell
$ echo $PATH
```
to check the list of your default PATH

3. Prepare bowtie2 index 
We have provided a pre-built bowtie2 index at our webpage for a quick start. Downloading requires some time (you can build the bowtie2 index yourself if you are using a genome version other than hg19 or mm9, see the Manual section). 

4. Prepare the peak calling software 
we use macs14 for peak calling in DrSeq2 on Drop-ChIP data and single cell ATAC-seq data. You can download MACS14 from here.

```shell
$ tar -zxvf MACS-1.4.2-1.tar.gz
$ cd MACS-1.4.2-1 
$ sudo python setup.py install # for root user 
$ python setup.py install --prefix /home/ # if you are not a root user
$ export PYTHONPATH=/home/lib/python2.6/site-packages:$PYTHONPATH
$ export PATH=/home/bin:$PATH More detailed about MACS see MACS web page
```
Step3.Run DrSeq2 on single cell ATAC-seq data
-------------------------------------------------------------------------------------------------------------------------------------------
You can run DrSeq2 pipeline to generate QC and analysis reports of your single cell ATAC-seq(scATAC) datasets.
Here, we provide an example of our simple mode on combined published scATAC-seq datasets (GSM1596255-GSM1596350, GSM1596735-GSM1596830 and GSM1597119-GSM1597214) and display DrSeq2 output in the following panel.
```shell
$ ATAC -i DrSeq2_scATAC -o scATAC -g hs --layout pair --geneannotation /PATHtoRefGene/hg19.refGene --cell_cutoff 5 --peak_cutoff 5 -C 3
```
For a brief description of major parameters, see the Manual section for more information.
-i input mapped reads file of single cell ATAC-seq. Sam/Bam files are supported,multiple files should seperate by comma.This parameter is required. <br>
-o output name.This parameter is required.Default is 'out'. <br>
-g genome_type for effective genome size. It can be shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9) Default:hs <br>
--layout pair-end sequencing or single-end sequencing.'pair' strands for pair-end sequencing.'single' strands for single-end sequencing. <br>
--geneannotation gene annotation file for CEAS. <br>
--cell_cutoff discard the peaks containing cells less than cell_cutoff. <br>
--peak_cutoff discard the cells containing peaks less than peak_cutoff. <br>
-C Given cluster numbers. DEFAULT: 3. <br>