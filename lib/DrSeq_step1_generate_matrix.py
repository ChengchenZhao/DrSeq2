#!/usr/bin/env python2.7

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import time
# --------------------------
# custom package
# --------------------------

### tool function
from DrSeq2pipe.Utility      import (Get,
                                   ChangeName,
                                   RaiseError,
                                   Log,
                                   LogError,
                                   LogCommand,
                                   Run,
                                   CreateDirectory,
                                   SampleDownTransformSam,
                                   GeneAnnotation,
                                   ReformBarcodeFastq,
                                   CombineReads,
                                   GenerateMatrix
                                   )
# --------------------------
# main 
# --------------------------

def step1_generate_matrix(conf_dict,logfile):
    '''
    generate expression matrix file 
    main data processing step, including mapping, generate expression matrix and QC matrix which is used in next step
    for fastq format : 
        STAR/bowtie2 mapping
        q30 filter, 
    for sam format:
        q30 filter     
    ''' 
    Log("Step1: alignment",logfile)
    t= time.time()
    ### create mapping dir 
    mapping_dir = conf_dict['General']['outputdirectory'] + 'mapping/'
    CreateDirectory(mapping_dir)
    ### check reads file format , start mapping step if format is fastq
    if conf_dict['General']['format'] == 'sam':
        Log('reads file format is sam, skip mapping step',logfile)
        conf_dict['General']['sam'] = conf_dict['General']['reads_file']
    else:
        Log('Now start mapping in %s , all mapping result will be here'%(mapping_dir),logfile)
        os.chdir(mapping_dir)
        ## choose mapping tool from STAR and bowtie2 according to config file
        if conf_dict['Step1_Mapping']['mapping_software_main'] == "STAR":
            Log('user choose STAR as alignment software',logfile)
            if Get('which STAR')[0].strip() == "":
                LogError('STAR is not detected in default PATH, make sure you installed STAR and export it into default PATH',logfile)
            mapping_cmd = 'STAR --genomeDir %s --readFilesIn %s --runThreadN %s'%(conf_dict['Step1_Mapping']['mapindex'],conf_dict['General']['reads_file'],conf_dict['Step1_Mapping']['mapping_p'])
            mapping_cmd2 = 'mv Aligned.out.sam %s.sam'%(conf_dict['General']['outname'])
            LogCommand(mapping_cmd,logfile)
            LogCommand(mapping_cmd2,logfile)
            
        elif conf_dict['Step1_Mapping']['mapping_software_main'] == "bowtie2":
            Log('user choose bowtie2 as alignment software',logfile)
            if Get('which bowtie2')[0].strip() == "":
                LogError('bowtie2 is not detected in default PATH, make sure you installed bowtie2 and export it into default PATH',logfile)
            mapping_cmd = 'bowtie2 -p %s -x %s -U %s -S %s.sam   2>&1 >>/dev/null |tee -a %s.bowtieout'%(conf_dict['Step1_Mapping']['mapping_p'],conf_dict['Step1_Mapping']['mapindex'],conf_dict['General']['reads_file'],conf_dict['General']['outname'],conf_dict['General']['outname'])
            LogCommand(mapping_cmd,logfile)
        
        elif conf_dict["Step1_Mapping"]["mapping_software_main"] == "HISAT2":
            Log('user choose HISAT2 as alignment software',logfile)
            if Get('which hisat2')[0].strip() == "":
                LogError('hisat2 is not detected in default PATH, make sure you installed hisat2 and export it into default PATH',logfile)
            mapping_cmd = 'hisat2 -p %s -x %s -U %s -S %s.sam   2>&1 >>/dev/null |tee -a %s.hisat2out'%(conf_dict['Step1_Mapping']['mapping_p'],conf_dict['Step1_Mapping']['mapindex'],conf_dict['General']['reads_file'],conf_dict['General']['outname'],conf_dict['General']['outname'])
            LogCommand(mapping_cmd,logfile)
        
        else:
            LogError("alignment tools can only be HISAT2, STAR or bowtie2",logfile)

        conf_dict['General']['sam'] = mapping_dir + conf_dict['General']['outname'] + '.sam'
    ### transform to bed file, awk helps to conduct q30 filtering
    Log("transfer sam file to aligned bed file with own script",logfile)
    conf_dict['General']['bed'] = mapping_dir + conf_dict['General']['outname'] + '.bed' 
    conf_dict['General']['sampledownsam'] = mapping_dir + conf_dict['General']['outname'] + '_sampledown.sam' 
    conf_dict['General']['sampledownbed'] = mapping_dir + conf_dict['General']['outname'] + '_sampledown.bed' 
    if int(conf_dict['Step1_Mapping']['q30filter']) == 1:
        Log("q30 filter is turned on",logfile)
    else:
        Log("q30 filter is turned off",logfile)
    ### use own script to transform sam to bed, and random sampling 5M mappable reads
    SampleDownTransformSam(conf_dict['General']['sam'],conf_dict['General']['bed'],conf_dict['General']['sampledownsam'],conf_dict['General']['sampledownbed'],5000000,int(conf_dict['Step1_Mapping']['q30filter']))
#        q30cmd = """samtools view -q 30 -XS %s | awk '{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if (substr($2,1,1) == "r") print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr" && $5 > 30) {if ($2 == 16) print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr" && $5 > 30) {if ($2 == 16) print $3,$4-1,$4,$1,255,"-";else print $3,$4-1,$4,$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        LogCommand(q30cmd,logfile,conf_dict['General']['dryrun'])
#        q30cmd = """samtools view -XS %s | awk '{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if (substr($2,1,1) == "r") print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if ($2 == 16) print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if ($2 == 16) print $3,$4-1,$4+length($11),$1,255,"-";else print $3,$4-1,$4,$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        LogCommand(q30cmd,logfile,conf_dict['General']['dryrun'])
    if not os.path.isfile(conf_dict['General']['bed']) or os.path.getsize(conf_dict['General']['bed']) == 0:
        LogError('Alignment step / q30 filtering step failed, check your alignment parameter and samfile',logfile)
    s1time = time.time() -t
    Log("time for alignment: %s"%(s1time),logfile)
    Log("Step1: alignment DONE",logfile)

    ### create annotation dir and generate related annotation file
    t = time.time() 
    Log("Step2: transform expression matrix",logfile)
    Log('generate related annotation file with own script',logfile)
    annotation_dir = conf_dict['General']['outputdirectory'] + 'annotation/'
    CreateDirectory(annotation_dir)
    os.chdir(annotation_dir)    
    GeneAnnotation(conf_dict['General']['gene_annotation'],conf_dict['Step2_ExpMat']['ttsdistance'],conf_dict['General']['outname'])

    ### create expression matrix dir and generate matrix
    Log('generate expression matrix and individual cell qc matrix with own script',logfile)    
    expdir = conf_dict['General']['outputdirectory'] + 'expmatrix/'
    CreateDirectory(expdir)
    os.chdir(expdir)
    
    ### use bedtools(intersect function) to assign exon/intron/intergenic/overlapping gene  information to all reads
    ### sort according to name
    Log('add gene annotation on aligned bed file',logfile)
    cmd1 = "bedtools intersect -a %s -b %s  -wo | sort -k 4,4 --parallel=6 -T . -S 8%% - > %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_symbol.bed',conf_dict['General']['outname']+'_on_symbol.bed')
    cmd2 = "bedtools intersect -a %s -b %s -c | sort -k 4,4 --parallel=6 -T . -S 8%% - > %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_cds.bed',conf_dict['General']['outname']+'_on_cds.bed')
    cmd3 = "bedtools intersect -a %s -b %s -c | sort -k 4,4 --parallel=6 -T . -S 8%% - > %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_3utr.bed',conf_dict['General']['outname']+'_on_3utr.bed')
    cmd4 = "bedtools intersect -a %s -b %s -c | sort -k 4,4 --parallel=6 -T . -S 8%% - > %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_5utr.bed',conf_dict['General']['outname']+'_on_5utr.bed')
    cmd5 = "bedtools intersect -a %s -b %s -c | sort -k 4,4 --parallel=6 -T . -S 8%% - > %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_TTSdis.bed',conf_dict['General']['outname']+'_on_TTSdis.bed')
    LogCommand(cmd1,logfile)
    LogCommand(cmd2,logfile)
    LogCommand(cmd3,logfile)
    LogCommand(cmd4,logfile)
    LogCommand(cmd5,logfile)

    ### transform barcode fastq to 3column txt file [name,cell_barcode,umi]
    if conf_dict['General']['format1'] == 'txt':
        Log('barcode files is reformed txt format, skip reform step',logfile)
        conf_dict['General']['barcode_reform'] = conf_dict['General']['barcode_file']
    else:
        Log('reform barcode files with own script',logfile)
        conf_dict['General']['barcode_reform'] = expdir + conf_dict['General']['outname'] + '_barcode_reform.txt'
        ReformBarcodeFastq(conf_dict['General']['barcode_file'],conf_dict['General']['barcode_reform'],conf_dict['General']['cell_barcode_range'],conf_dict['General']['umi_range'])
    ### sort according name
    cmdsort = 'sort -k 1,1 --parallel=6 -T . -S 8%% %s > %s'%(conf_dict['General']['barcode_reform'],expdir + conf_dict['General']['outname'] + '_barcode_reform_sort.txt')
    LogCommand(cmdsort,logfile)
    conf_dict['General']['barcode_reform'] = expdir + conf_dict['General']['outname'] + '_barcode_reform_sort.txt'
    
    ### combine gene annotation, reads, barcode together
    Log('combine annotation and barcode on reads with own script',logfile)
    CombineReads(conf_dict['General']['barcode_reform'],conf_dict['General']['outname']+'_on_cds.bed',conf_dict['General']['outname']+'_on_3utr.bed',conf_dict['General']['outname']+'_on_5utr.bed',conf_dict['General']['outname']+'_on_symbol.bed',conf_dict['General']['outname']+'_on_TTSdis.bed',conf_dict['General']['outname']+ '_combined.bed',conf_dict['Step2_ExpMat']['duplicate_measure'])   
     
    ### sort combined file by umi+loci, for following duplicate detection
    cmd6 = "sort -k 7,7 -k 5,5 --parallel=6 -T . -S 8%% %s > %s"%(conf_dict['General']['outname']+ '_combined.bed',conf_dict['General']['outname']+ '_combined_sort.bed')
    LogCommand(cmd6,logfile)
    
    ### generate expression and QC matrix based on combined file
    Log('generate expression matrix and QC matrix with own script',logfile)
    ### qcmatfull contains all cell_barcodes, while qcmat,expmat only contain cell_barcodes >= covergncutoff(100, default)
    conf_dict['Step2_ExpMat']['qcmatfull'] = expdir + conf_dict['General']['outname'] + "_qcmatfull.txt"    
    conf_dict['Step2_ExpMat']['qcmat'] = expdir + conf_dict['General']['outname'] + "_qcmat.txt"
    conf_dict['Step2_ExpMat']['expmat'] = expdir + conf_dict['General']['outname'] + "_expmat.txt"
    
    GenerateMatrix(conf_dict['General']['gene_annotation'],conf_dict['General']['outname']+ '_combined_sort.bed',conf_dict['Step2_ExpMat']['filterttsdistance'],conf_dict['Step2_ExpMat']['qcmatfull'],conf_dict['Step2_ExpMat']['qcmat'],conf_dict['Step2_ExpMat']['expmat'],conf_dict['Step2_ExpMat']['covergncutoff'],conf_dict['Step2_ExpMat']['umidis1'])    
        
    Log("Step2 transform expression matrix DONE",logfile)
    s2time = time.time() -t
    Log("time for transform expmat: %s"%(s2time),logfile)
    conf_dict['results'] = {}
    #conf_dict['results']['expmat'] = conf_dict['Step2_ExpMat']['expmat']
    #conf_dict['results']['qcmat'] = conf_dict['Step2_ExpMat']['qcmat']
    
    return conf_dict
