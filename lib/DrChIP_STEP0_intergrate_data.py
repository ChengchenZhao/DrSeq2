#!/usr/bin/env python

# ================================================================================
# Python Modual
# ================================================================================

import os
import sys
from string import *

# ================================================================================
# custom package
# ================================================================================


### tool function
from DrSeq2pipe.Utility import (Get,GetError,DetectMemory,LogError,Log)

# ================================================================================
# function 
# ================================================================================

# input : barcode list , fastq_1 ,fastq_2

def Step0IntegrateData(conf_dict,logfile):
    '''
    step0 integrate data 
    check and complement parameter
    '''
    Log("Start Dr.ChIP",logfile)
    Log("Step0: Data integrate",logfile)
    
    ### check output name
    if "/" in conf_dict['General']['outname']:
        LogError("outname is the name of all your output result, cannot contain "/", current outname is  %s"%(conf_dict['General']['outname']),logfile)
    ### check data path , format ,
    if "~" in conf_dict['General']['fastq_1']:
        conf_dict['General']['fastq_1'] = os.path.expanduser(conf_dict['General']['fastq_1'])
    if "~" in conf_dict['General']['fastq_2']:
        conf_dict['General']['fastq_2'] = os.path.expanduser(conf_dict['General']['fastq_2'])
    if "~" in conf_dict['General']['barcode_file']:
        conf_dict['General']['barcode_file'] = os.path.expanduser(conf_dict['General']['barcode_file'])
    if not conf_dict['General']['fastq_1'].startswith('/'):
        conf_dict['General']['fastq_1'] = conf_dict['General']['startdir'] + conf_dict['General']['fastq_1']
    if not conf_dict['General']['fastq_2'].startswith('/'):
        conf_dict['General']['fastq_2'] = conf_dict['General']['startdir'] + conf_dict['General']['fastq_2']
    if not conf_dict['General']['barcode_file'].startswith('/'):
        conf_dict['General']['barcode_file'] = conf_dict['General']['startdir'] + conf_dict['General']['barcode_file']
    if not os.path.isfile(conf_dict['General']['fastq_1']):
        LogError("fastq_1 file %s not found"%(conf_dict['General']['fastq_1']),logfile)
    if not os.path.isfile(conf_dict['General']['fastq_2']):
        LogError("fastq_2 file %s not found"%(conf_dict['General']['fastq_2']),logfile)
    if not os.path.isfile(conf_dict['General']['barcode_file']):
        LogError("barcode_file file %s not found"%(conf_dict['General']['barcode_file']),logfile)
        
    if not (conf_dict['General']['fastq_1'].endswith('.fastq') and conf_dict['General']['fastq_2'].endswith('.fastq')):
        LogError("input files should be fastq files.",logfile)
    else:
        Log('Detected input file format is fastq',logfile)
        conf_dict['General']['format'] = 'fastq'

    ### check barcode range
    conf_dict['General']['barcode_range_1'] = conf_dict['General']['barcode_range_1']
    conf_dict['General']['barcode_range_2'] = conf_dict['General']['barcode_range_2']
    conf_dict['General']['barcode_file_range_1'] = conf_dict['General']['barcode_file_range_1']
    conf_dict['General']['barcode_file_range_2'] = conf_dict['General']['barcode_file_range_2']
    ### check gene annotation file
    if conf_dict['General']['gene_annotation'] == "":
        LogError("gene annotation file cannot be empty",logfile)
    if not "/" in conf_dict['General']['gene_annotation'] : 
        LogError("absolute path for gene annotation file required",logfile)        
    if not os.path.isfile(conf_dict['General']['gene_annotation'] ):
        LogError("cannot find gene annotation file : %s"%(conf_dict['General']['gene_annotation'] ),logfile)
        
    ### mapping index
    conf_dict['Step1_Mapping']['mapindex'] = os.path.expanduser(conf_dict['Step1_Mapping']['mapindex'])
    if conf_dict['General']['format'] == 'fastq':
        if conf_dict['Step1_Mapping']['mapping_software'] == "bowtie2":
            Log('use bowtie2 as alignment tools',logfile)
#            conf_dict['Step1_Mapping']['mapindex'] = indexdir + conf_dict['General']['genome_version']
            indexfile1 = conf_dict['Step1_Mapping']['mapindex']+'.1.bt2'
            if not os.path.isfile(indexfile1):
                LogError("cannot find bowtie2 index file : %s "%(indexfile1),logfile)
        else:
            LogError("alignment tools can only be bowtie2 by now",logfile)


    ### check options
    Log('option setting: ',logfile)
    try:
        Log('mapping thread is %s'%(str(int(conf_dict['Step1_Mapping']['p']))),logfile)
    except:
        LogError('p should be int, current value is %s'%(conf_dict['Step1_Mapping']['p']),logfile)
        
    if not int(conf_dict['Step1_Mapping']['q30filter']) in [0,1]:
        LogError('q30filter measurement can only be 0/1, current value is %s'%(conf_dict['Step1_Mapping']['q30filter']),logfile)


    ### check Rscript
    if not 'Usage' in GetError('Rscript')[1] and not 'version' in GetError('Rscript')[1]:
       LogError('require Rscript',logfile)
    
    ### check pdflatex
    if Get('pdflatex --help')[0] == "":
        Log('pdflatex was not installed, Dr.ChIP is still processing but no summary QC report generated',logfile)
        conf_dict['General']['latex'] = 0
    else:
        conf_dict['General']['latex'] = 1

    Log('Step0 Data integrate DONE',logfile)

    return conf_dict
