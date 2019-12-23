#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
from string import *

# --------------------------
# custom package
# --------------------------


### tool function
import DrSeq2pipe
from DrSeq2pipe.Utility      import *
# --------------------------
# main 
# --------------------------

def step0_integrate_data(conf_dict,logfile):
    '''
    step0 integrate data 
    check and complement parameter
    '''
    Log("Start Drseq",logfile)
    Log("Step0: Data integrate",logfile)
    
    ### check output name
    if "/" in conf_dict['General']['outname']:
        LogError("outname is the name of all your output result, cannot contain "/", current outname is  %s"%(conf_dict['General']['outname']),logfile)
    ### check data path , format ,
    if "~" in conf_dict['General']['barcode_file']:
        conf_dict['General']['barcode_file'] = os.path.expanduser(conf_dict['General']['barcode_file'])
    if "~" in conf_dict['General']['reads_file']:
        conf_dict['General']['barcode_file'] = os.path.expanduser(conf_dict['General']['barcode_file'])
    if not conf_dict['General']['barcode_file'].startswith('/'):
        conf_dict['General']['barcode_file'] = conf_dict['General']['startdir'] + conf_dict['General']['barcode_file']
    if not conf_dict['General']['reads_file'].startswith('/'):
        conf_dict['General']['reads_file'] = conf_dict['General']['startdir'] + conf_dict['General']['reads_file']
    
    if not os.path.isfile(conf_dict['General']['barcode_file']):
        LogError("barcode file %s not found"%(conf_dict['General']['barcode_file']),logfile)
    if not os.path.isfile(conf_dict['General']['reads_file']):
        LogError("reads file %s not found"%(conf_dict['General']['reads_file']),logfile)
        
    if not conf_dict['General']['barcode_file'].endswith('.fastq'):
        if conf_dict['General']['barcode_file'].endswith('.txt'):
            Log('barcode file is reformed txt file',logfile)
            conf_dict['General']['format1'] = 'txt'   
        else:
            LogError("barcode file is not a fastq file: %s"%(conf_dict['General']['barcode_file']),logfile)
    else:
        conf_dict['General']['format1'] = 'fastq'
    if conf_dict['General']['reads_file'].endswith('.fastq') or conf_dict['General']['reads_file'].endswith('.fq'):
        conf_dict['General']['format'] = 'fastq'
        Log('Detected input file format is fastq',logfile)
    elif conf_dict['General']['reads_file'].endswith('.sam'): 
        conf_dict['General']['format'] = 'sam'
        Log('Detected input file format is sam',logfile)
    else:
        LogError("reads file is not a fastq or sam file: %s"%(conf_dict['General']['reads_file']),logfile)
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
#        if not conf_dict['Step1_Mapping']['mapindex'].endswith("/"):
#            conf_dict['Step1_Mapping']['mapindex'] += "/"
        if conf_dict['Step1_Mapping']['mapping_software_main'] == "STAR":
            Log('use STAR as alignment tools',logfile)
            if int(conf_dict['Step1_Mapping']['checkmem']) == 1:
                Log('memory check is turned on, check total memory',logfile)
                totalMemory = DetectMemory()
                if totalMemory == "NA"  : 
                    LogError('''cannot detect total memory (because your server don't have /proc/meminfo file or you are running Dr.seq on Mac computer), Dr.seq exit to protect your server from crash down. You can turn off the memory check and run Dr.seq again if you do want to use STAR as mapping software or you can use bowtie2 instead.''',logfile)
                elif totalMemory < 40 : 
                    LogError('''Total memory of your server/computer is %sG, less than 40G (memory cutoff for STAR), Dr.seq exit to protect your server from crash down. You can turn off the memory check and run Dr.seq again if you do want to use STAR as mapping software or you can use bowtie2 instead '''%(str(totalMemory)),logfile)
                else:
                    Log('''Total memory of your  server/computer is %sG, greater than 40G (memory cutoff for STAR), Dr.seq will use STAR as mapping software'''%(str(totalMemory)),logfile)
            else:
                Log('memory check is turned off, start mapping with STAR ### STAR consume > 30G memory, make sure your server have enough memory ###',logfile)
#            conf_dict['Step1_Mapping']['mapindex'] +='%s.star'%(conf_dict['General']['genome_version'])
            if not os.path.isdir(conf_dict['Step1_Mapping']['mapindex']):
                LogError("cannot find STAR index folder : %s"%(conf_dict['Step1_Mapping']['mapindex']),logfile)
        elif conf_dict['Step1_Mapping']['mapping_software_main'] == "bowtie2":
            Log('use bowtie2 as alignment tools',logfile)
#            conf_dict['Step1_Mapping']['mapindex'] = indexdir + conf_dict['General']['genome_version']
            indexfile1 = conf_dict['Step1_Mapping']['mapindex']+'.1.bt2'
#           if not os.path.isdir(indexdir):
#               LogError("cannot find bowtie2 index folder : %s "%(indexdir),logfile)
            if not os.path.isfile(indexfile1):
                LogError("cannot find bowtie2 index file : %s "%(indexfile1),logfile)
        elif conf_dict['Step1_Mapping']['mapping_software_main'] == "HISAT2":
            Log('use HISAT2 as alignment tools',logfile)
#            conf_dict['Step1_Mapping']['mapindex'] = indexdir + conf_dict['General']['genome_version']
            indexfile1 = conf_dict['Step1_Mapping']['mapindex']+'.1.ht2'
#           if not os.path.isdir(indexdir):
#               LogError("cannot find bowtie2 index folder : %s "%(indexdir),logfile)
            if not os.path.isfile(indexfile1):
                LogError("cannot find HISAT2 index file : %s "%(indexfile1),logfile)
        else:
            LogError("alignment tools can only be HISAT2, STAR and bowtie2",logfile)


    ### check options
    Log('option setting: ',logfile)
    try:
        Log('mapping thread is %s'%(str(int(conf_dict['Step1_Mapping']['mapping_p']))),logfile)
    except:
        LogError('mapping_p should be int, current value is %s'%(conf_dict['Step1_Mapping']['mapping_p']),logfile)
        
    if not int(conf_dict['Step1_Mapping']['q30filter']) in [0,1]:
        LogError('q30filter measurement can only be 0/1, current value is %s'%(conf_dict['Step1_Mapping']['q30filter']),logfile)

    if not int(conf_dict['Step2_ExpMat']['filterttsdistance']) in [0,1]:
        LogError('filterttsdistance measurement can only be 0/1, current value is %s'%(conf_dict['Step2_ExpMat']['filterttsdistance']),logfile)
    
    if not int(conf_dict['Step2_ExpMat']['ttsdistance']) > 0:
        LogError('ttsdistance value should greater than 0, current value is %s'%(conf_dict['Step2_ExpMat']['ttsdistance']),logfile)

    if	 int(conf_dict['Step2_ExpMat']['covergncutoff']) > 10000:
        LogError('covergncutoff value cannot be greater than 10000, current value is %s'%(conf_dict['Step2_ExpMat']['covergncutoff']),logfile)
    
    if not int(conf_dict['Step2_ExpMat']['duplicate_measure']) in [0,1,2,3]:
        LogError('duplicate_measure value can only be 0~3, current value is %s'%(conf_dict['Step2_ExpMat']['duplicate_measure']),logfile)

    if not int(conf_dict['Step3_QC']['select_cell_measure']) in [1,2]:
        LogError('select_cell_measure value can only be 1 or 2, current value is %s'%(conf_dict['Step3_QC']['select_cell_measure']),logfile)
    
    if int(conf_dict['Step3_QC']['select_cell_measure']) == 1:
        try:
            int(conf_dict['Step3_QC']['covergncluster'])
        except:
            LogError('covergncluster value should be integer, current value is %s'%(conf_dict['Step3_QC']['covergncluster']),logfile)
    elif int(conf_dict['Step3_QC']['select_cell_measure']) == 2:
        try:
            int(conf_dict['Step3_QC']['topumicellnumber'])
        except:
            LogError('topumicellnumber value should be integer, current value is %s'%(conf_dict['Step3_QC']['covergncluster']),logfile)   
    else: 
        LogError('select_cell_measure value can only be 1 or 2, current value is %s'%(conf_dict['Step3_QC']['select_cell_measure']),logfile)

    if not int(conf_dict['Step3_QC']['remove_low_dup_cell']) in [0,1]:
        LogError('remove_low_dup_cell measurement can only be 0/1, current value is %s'%(conf_dict['Step3_QC']['remove_low_dup_cell']),logfile)
    if float(conf_dict['Step3_QC']['non_dup_cutoff']) <= 0  or float(conf_dict['Step3_QC']['non_dup_cutoff']) >=1 :
        LogError('non_dup_cutoff measurement should be in 0~1, current value is %s'%(conf_dict['Step3_QC']['non_dup_cutoff']),logfile)
    if float(conf_dict['Step4_Analysis']['highvarz']) <= 0  :
        LogError('non_dup_cutoff measurement cannot be <= 0, current value is %s'%(conf_dict['Step4_Analysis']['highvarz']),logfile)
    if float(conf_dict['Step4_Analysis']['selectpccumvar']) <= 0 or float(conf_dict['Step4_Analysis']['selectpccumvar']) >=1 :
        LogError('selectpccumvar measurement should be in 0~1, current value is %s'%(conf_dict['Step4_Analysis']['selectpccumvar']),logfile)
    if not int(conf_dict['Step4_Analysis']['clustering_method']) in [1,2,3,4] :
        LogError('clustering_method measurement should be chosen from 1,2,3 and 4, current value is %s'%(conf_dict['Step4_Analysis']['clustering_method']),logfile)
    if not int(conf_dict['Step4_Analysis']['dimensionreduction_method']) in [1,2,3] :
        LogError('dimensionreduction_method measurement should be chosen from 1,2 and 3 , current value is %s'%(conf_dict['Step4_Analysis']['dimensionreduction_method']),logfile)

    if int(conf_dict['Step4_Analysis']['dimensionreduction_method'])==3 and not os.path.isfile(os.path.join(DrSeq2pipe.__path__[0], "Rscript/projsplx_R.so")):
        LogError(""" 2 R packages : "Matrix", "parallel" and an external C code : "projsplx_R.so" are necessary for your selected dimension reduction method (SIMLR), please install these two R packages first. We provide supporting information in the manual file and our web page. Or you can choose other method for dimension reduction. """ ,logfile)


    ### check Rscript
    #if not 'Usage' in GetError('Rscript')[1] and not 'version' in GetError('Rscript')[1]:
    #    LogError('require Rscript',logfile)
    
    ### check pdflatex
    if Get('pdflatex --help')[0] == "":
        Log('pdflatex was not installed, Dr.seq is still processing but no summary QC report generated',logfile)
        conf_dict['General']['latex'] = 0
    else:
        conf_dict['General']['latex'] = 1

    Log('Step0 Data integrate DONE',logfile)



    return conf_dict
    
    
    
