#!/usr/bin/env python

# ================================================================================
# Python Modual
# ================================================================================

import os
import sys
import string
import time

# ================================================================================
# custom package
# ================================================================================


### tool function
from DrSeq2pipe.Utility import *

# ================================================================================
# function 
# ================================================================================

def ComplementationPairingSeq(sequence):
    complementation_pairing_dic = {"A":"T","C":"G","T":"A","G":"C","N":"N"}
    sequence_pair = ""
    for i in sequence[::-1]:
        sequence_pair += complementation_pairing_dic[i]
    return sequence_pair

def BarcodeDictionary(barcode_file,barcode_file_range_1,barcode_file_range_2):
    barcode_info = open(barcode_file)
    barcode_dic_1 = {}
    barcode_dic_2 = {}
    for each in barcode_info:
        each = each.strip().split()
        cell_id = each[0]
        barcode_1 = each[1][(int(barcode_file_range_1.split(":")[0])-1):int(barcode_file_range_1.split(":")[1])]
        barcode_2 = each[1][(int(barcode_file_range_2.split(":")[0])-1):int(barcode_file_range_2.split(":")[1])]
        barcode_dic_1[barcode_1] = cell_id
        barcode_dic_1[ComplementationPairingSeq(barcode_1)] = cell_id
        barcode_dic_2[barcode_2] = cell_id
        barcode_dic_2[ComplementationPairingSeq(barcode_2)] = cell_id
    barcode_info.close()
    return barcode_dic_1,barcode_dic_2

def BarcodeFilteredReads(barcode_dic_1,barcode_dic_2,fastq_file1,fastq_file2,barcode_range_1,barcode_range_2,logfile,mapping_dir,outname,conf_dict):
    fastq1_info = open(fastq_file1)
    fastq2_info = open(fastq_file2)
    BC1_reads_dic = {}
    BC2_reads_dic = {}
    BC_filter_reads_dic = {}

    count=1
    for each in fastq1_info:
        if count%4 == 1:
            read_id = each.strip().split(" ")[0][1:]
        if count%4 == 2:
            bc_seq = each.strip()[(int(barcode_range_1.split(":")[0])-1):int(barcode_range_1.split(":")[1])]
            if barcode_dic_1.has_key(bc_seq):
                BC1_reads_dic[read_id] = barcode_dic_1[bc_seq]
        count+=1
    Log("number of all reads :\t%d\tnumber of BC1 reads :\t%d"%(count/4,len(BC1_reads_dic)),logfile)

    count=1
    for each in fastq2_info:
        if count%4 == 1:
            read_id = each.strip().split(" ")[0][1:]
        if count%4 == 2:
            bc_seq = each.strip()[(int(barcode_range_2.split(":")[0])-1):int(barcode_range_2.split(":")[1])]
            if barcode_dic_2.has_key(bc_seq):
                BC2_reads_dic[read_id] = barcode_dic_2[bc_seq]
        count+=1
    Log("number of all reads :\t%d\tnumber of BC2 reads :\t%d"%(count/4,len(BC2_reads_dic)),logfile)

    for each_read in set(BC1_reads_dic.keys())&set(BC2_reads_dic.keys()):
        if BC1_reads_dic[each_read] == BC2_reads_dic[each_read]:
            BC_filter_reads_dic[each_read] = BC1_reads_dic[each_read]
    
    Log("number of reads after BC filtered:\t%d"%(len(BC_filter_reads_dic)),logfile)
    fastq1_info.close()
    fastq2_info.close()
    
    conf_dict['Step2_QC']['total_reads_pair_N'] = count/4
    conf_dict['Step2_QC']['bc_filter_reads_pair_N'] = len(BC_filter_reads_dic)


    fastq1_info = open(fastq_file1)
    fastq2_info = open(fastq_file2)
    new_fastq1_info = open(mapping_dir + outname + "selected_1.fastq","w")
    new_fastq2_info = open(mapping_dir + outname + "selected_2.fastq","w")
    count=1
    for each in fastq1_info:
        if count%4 == 1:
            read_id = each.strip().split(" ")[0][1:]
            if BC_filter_reads_dic.has_key(read_id):
                selected = True
                line1 = each.strip().split(" ")[0]+":BARCODE"+str(BC_filter_reads_dic[read_id])+" "+" ".join(each.strip().split(" ")[1:])+"\n"
                # print line1
                new_fastq1_info.write(line1)
            else:
                selected = False
        else:
            if selected == True:
                new_fastq1_info.write(each)
        count+=1
    count=1
    for each in fastq2_info:
        if count%4 == 1:
            read_id = each.strip().split(" ")[0][1:]
            if BC_filter_reads_dic.has_key(read_id):
                selected = True
                line1 = each.strip().split(" ")[0]+":BARCODE"+str(BC_filter_reads_dic[read_id])+" "+" ".join(each.strip().split(" ")[1:])+"\n"
                new_fastq2_info.write(line1)
            else:
                selected = False
        else:
            if selected == True:
                new_fastq2_info.write(each)
        count+=1
    new_fastq1_info.close()
    new_fastq2_info.close()

    conf_dict['General']['fastq_1'] = mapping_dir + outname + "selected_1.fastq"
    conf_dict['General']['fastq_2'] = mapping_dir + outname + "selected_2.fastq"

    return conf_dict,BC_filter_reads_dic,BC1_reads_dic

def lengthFilterReads(sam_file,BC1_reads_dic,outname,conf_dict,out_dir,logfile):
    
    filter_length = 1000

    command = "samtools view -SH %s.sam > %s_head.sam"%(outname,outname)
    LogCommand(command,logfile)
    command = "samtools view -S %s.sam | awk '{if ($9 <= 1000 && $9 >-1000 && $9!=0) print $0}' > %s_selected_reads.sam"%(outname,outname)
    LogCommand(command,logfile)
    command = "cat %s_head.sam %s_selected_reads.sam > %s.selected.sam"%(outname,outname,outname)
    LogCommand(command,logfile)
    command = "samtools view -Sb %s.selected.sam > %s.selected.bam"%(outname,outname)
    LogCommand(command,logfile)
    command = """bamToBed -bedpe -i %s.selected.bam | awk '{if($2<$5) {split($7, array, "BARCODE");print $1"\t"$2"\t"$6"\t"array[2];} if($2>$5) {split($7, array, "BARCODE");print $1"\t"$5"\t"$3"\t"array[2]}}' | sort -k1,1 -k2,2n > %s.selected.bed"""%(outname,outname)
    LogCommand(command,logfile)
    
    mapped_paireds = int(os.popen("samtools flagstat %s.bam"%(conf_dict['General']['outputdirectory'] + 'mapping/' +conf_dict['General']['outname'])).read().split("read2")[1].split("+")[0])/2
    len_filtered_reads_number = os.popen("wc -l %s.selected.bed"%outname).read().split("%s.selected.bed"%outname)[0].replace(" ","")
    Log("Number of all mapped reads pair:\t%d" % mapped_paireds , logfile)
    Log("Number of final reads pair after length filter :\t%s"%len_filtered_reads_number,logfile)
    
    conf_dict['Step2_QC']['mapped_reads_pair_N'] = mapped_paireds
    conf_dict['Step2_QC']['final_reads_pair_N'] = len_filtered_reads_number

    command = "cut -f 4 %s.selected.bed | sort | uniq -c > %s.readsnumber.txt"%(outname,outname)
    LogCommand(command,logfile)
    conf_dict['General']['sam'] = out_dir+"%s.selected.sam"%outname
    return conf_dict

def Step1Mapping(conf_dict,logfile):
    '''
    generate expression matrix file 
    main data processing step, including barcode deconvolution, mapping, and QC matrix which is used in next step
    for fastq format : 
        STAR/bowtie2 mapping
        q30 filter
    ''' 
    Log("[1] barcode deconvolution",logfile)
    t = time.time()
    ### create selected fastq dir
    mapping_dir = conf_dict['General']['outputdirectory'] + 'fastq/'
    CreateDirectory(mapping_dir)
    barcode_dic_1,barcode_dic_2 = BarcodeDictionary(conf_dict['General']['barcode_file'],conf_dict['General']['barcode_file_range_1'],conf_dict['General']['barcode_file_range_2'])
    conf_dict,BC_filter_reads_dic,BC1_reads_dic = BarcodeFilteredReads(barcode_dic_1,barcode_dic_2,conf_dict['General']['fastq_1'],conf_dict['General']['fastq_2'],conf_dict['General']['barcode_range_1'],conf_dict['General']['barcode_range_2'],logfile,mapping_dir,conf_dict['General']['outname'],conf_dict)

    Log("[2] alignment",logfile)
    t= time.time()
    ### create mapping dir 
    mapping_dir = conf_dict['General']['outputdirectory'] + 'mapping/'
    CreateDirectory(mapping_dir)
    ### check reads file format , start mapping step if format is fastq
    if conf_dict['General']['format'] == 'sam':
        Log('reads file format is sam, skip mapping step',logfile)
        conf_dict['General']['sam'] = conf_dict['General']['sam_file']
    else:
        Log('Now start mapping in %s , all mapping result will be here'%(mapping_dir),logfile)
        os.chdir(mapping_dir)
        ## choose mapping tool from STAR and bowtie2 according to config file
        if conf_dict['Step1_Mapping']['mapping_software'] == "bowtie2":
            Log('user choose bowtie2 as alignment software',logfile)
            if Get('which bowtie2')[0].strip() == "":
                LogError('bowtie2 is not detected in default PATH, make sure you installed bowtie2 and export it into default PATH',logfile)
            mapping_cmd = 'bowtie2 -X %s --trim5 %s -p %s -x %s -1 %s -2 %s -S %s.sam 2>&1 >>/dev/null |tee -a %s.bowtieout'%(conf_dict['Step1_Mapping']['x'],conf_dict['Step1_Mapping']['trim5'],conf_dict['Step1_Mapping']['p'],conf_dict['Step1_Mapping']['mapindex'],conf_dict['General']['fastq_1'],conf_dict['General']['fastq_2'],conf_dict['General']['outname'],conf_dict['General']['outname'])
            LogCommand(mapping_cmd,logfile)
            sam2bam_cmd = 'samtools view -Sb %s.sam > %s.bam'%(conf_dict['General']['outname'],conf_dict['General']['outname'])
            LogCommand(sam2bam_cmd,logfile)
        # elif conf_dict['Step1_Mapping']['mapping_software'] == "STAR":
        #     Log('user choose STAR as alignment software',logfile)
        #     if Get('which STAR')[0].strip() == "":
        #         LogError('STAR is not detected in default PATH, make sure you installed STAR and export it into default PATH',logfile)
        #     mapping_cmd = 'STAR --genomeDir %s --readFilesIn %s --runThreadN %s'%(conf_dict['Step1_Mapping']['mapindex'],conf_dict['General']['fastq_1'],conf_dict['Step1_Mapping']['p'])
        #     mapping_cmd2 = 'mv Aligned.out.sam %s.sam'%(conf_dict['General']['outname'])
        #     LogCommand(mapping_cmd,logfile)
        #     LogCommand(mapping_cmd2,logfile)
        else:
            LogError("alignment tools can only be bowtie2",logfile)

        conf_dict['General']['sam'] = mapping_dir + conf_dict['General']['outname'] + '.sam'
    
    ## filter the reads that more than 1000 between two pair reads.
    conf_dict = lengthFilterReads(conf_dict['General']['sam'],BC1_reads_dic,conf_dict['General']['outname'],conf_dict,mapping_dir,logfile)

    ### transform to bed file, awk helps to conduct q30 filtering
    Log("transfer sam file to aligned bed file with own script",logfile)
    conf_dict['General']['bed'] = mapping_dir + conf_dict['General']['outname'] + '.bed' 
    conf_dict['General']['sampledownsam'] = mapping_dir + conf_dict['General']['outname'] + '_sampledown.sam' 
    conf_dict['General']['sampledownbed'] = mapping_dir + conf_dict['General']['outname'] + '_sampledown.bed' 
    if int(conf_dict['Step1_Mapping']['q30filter']) == 1:
        Log("q30 filter is turned on",logfile)
        ## use own script to transform sam to bed, and random sampling 5M mappable reads
        SampleDownTransformSam(conf_dict['General']['sam'],conf_dict['General']['bed'],conf_dict['General']['sampledownsam'],conf_dict['General']['sampledownbed'],5000000,int(conf_dict['Step1_Mapping']['q30filter']))
        q30cmd = """samtools view -q 30 -XS %s | awk '{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if (substr($2,1,1) == "r") print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr" && $5 > 30) {if ($2 == 16) print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr" && $5 > 30) {if ($2 == 16) print $3,$4-1,$4,$1,255,"-";else print $3,$4-1,$4,$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
        LogCommand(q30cmd,logfile)
    else:
        Log("q30 filter is turned off",logfile)
        q30cmd = """samtools view -XS %s | awk '{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if (substr($2,1,1) == "r") print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if ($2 == 16) print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if ($2 == 16) print $3,$4-1,$4+length($11),$1,255,"-";else print $3,$4-1,$4,$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
        LogCommand(q30cmd,logfile)
    if not os.path.isfile(conf_dict['General']['bed']) or os.path.getsize(conf_dict['General']['bed']) == 0:
        LogError('Alignment step / q30 filtering step failed, check your alignment parameter and samfile',logfile)
    s1time = time.time() - t
    Log("time for alignment: %s"%(s1time),logfile)
    Log("Step1: alignment DONE",logfile)
    
    return conf_dict

    