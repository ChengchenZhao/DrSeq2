#!/usr/bin/env python2.7

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

def Step2_QC(conf_dict,logfile):
    '''
    start RseQC
    mapping stat
    single cell level QC
    '''
    # start
    # create section for 
    
    Log('Step2: bulk and individual cell QC',logfile)
    
    ## reads quality
    t= time.time()
    
    QC_dir = conf_dict['General']['outputdirectory'] + 'QC/'
    CreateDirectory(QC_dir)
    os.chdir(QC_dir)

    readsqc(conf_dict['General']['sampledownsam'],conf_dict['General']['outname'])
    Log('generate bulk cell QC measurement with own script, based on sample down reads',logfile)

    cmd = "%s %s %s"%('Rscript',conf_dict['rscript']+'DrChIP_readsbulk_QC.r',conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    ## peak calling based on bulk cell
    Log('call peaking based on combined reads using macs14',logfile)
    if Get('which macs14')[0].strip() == "":
        LogError('macs14 is not detected in default PATH, make sure you installed macs14 and export it into default PATH',logfile)
    # cmd = "macs14 -g %s -SB -p %s --keep-dup 1 --nomodel --shiftsize 73 -t %s -n %s"%(conf_dict['General']["species"],conf_dict["Step2_QC"]["peak_calling_cutoff"],conf_dict['General']['sam'],conf_dict['General']['outname'])
    # LogCommand(cmd,logfile)
    # cmd = "gunzip %s_MACS_bedGraph/treat/%s_treat_afterfiting_all.bdg.gz"%(conf_dict['General']['outname'],conf_dict['General']['outname'])
    # LogCommand(cmd,logfile)
    # cmd = "sort -k1,1 -k2,2n %s_MACS_bedGraph/treat/%s_treat_afterfiting_all.bdg > %s.sort.bdg"%(conf_dict['General']['outname'],conf_dict['General']['outname'],conf_dict['General']['outname'])
    # LogCommand(cmd,logfile)
    # # time.sleep(200)
    # cmd = "bedClip %s.sort.bdg %s %s.tmp.bdg"%(conf_dict['General']['outname'],conf_dict['General']['chrom_length'],conf_dict['General']['outname'])
    # LogCommand(cmd,logfile)
    # cmd = "grep -v _ %s.tmp.bdg > %s.clean.bdg"%(conf_dict['General']['outname'],conf_dict['General']['outname'])
    # LogCommand(cmd,logfile)
    # cmd = "bedGraphToBigWig %s.clean.bdg %s %s.bw"%(conf_dict['General']['outname'],conf_dict['General']['chrom_length'],conf_dict['General']['outname'])
    # LogCommand(cmd,logfile)
    # cmd = "ceasBW -w %s.bw -b %s_peaks.bed -g %s -l %s --name %s_ceas > ceasBW.out"%(conf_dict['General']["outname"],conf_dict['General']["outname"],conf_dict['General']['gene_annotation'],conf_dict['General']['chrom_length'],conf_dict['General']['outname'])
    # LogCommand(cmd,logfile)
    cmd = """macs14 -g %s -Sw -p %s --keep-dup all --nomodel --shiftsize 73 -t %s -n %s"""%(conf_dict['General']["species"],conf_dict["Step2_QC"]["peak_calling_cutoff"],conf_dict['General']['sam'],conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    cmd = """gunzip %s_MACS_wiggle/treat/%s_treat_afterfiting_all.wig.gz""" %(conf_dict['General']['outname'],conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    cmd = """ceas -w %s_MACS_wiggle/treat/%s_treat_afterfiting_all.wig -b %s_peaks.bed -g %s --name %s_ceas""" %(conf_dict['General']['outname'],conf_dict['General']['outname'],conf_dict['General']['outname'],conf_dict['General']['gene_annotation'],conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    
    ## mapping ratio QC
    mapping_dir = conf_dict['General']['outputdirectory'] + 'mapping/'
    mapping_ratio = os.popen("samtools flagstat %s.bam"%(mapping_dir+conf_dict['General']['outname'])).read().split("properly paired ")[1].split("%")[0][1:]
    cmd = "Rscript %s %s %s"%(conf_dict['rscript']+'DrChIP_mappingRatio_cdf.r',mapping_ratio,conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    
    if not conf_dict.has_key('QCplots'):
        conf_dict['QCplots'] = {}
    if not conf_dict.has_key('results'):
        conf_dict['results'] = {}
    conf_dict['QCplots']['read_qul'] = QC_dir + conf_dict['General']['outname'] + '_Figure1_quality_heatmap.pdf'
    conf_dict['QCplots']['read_nvc'] = QC_dir + conf_dict['General']['outname'] + '_Figure2_NVC.pdf'
    conf_dict['QCplots']['read_gc'] = QC_dir + conf_dict['General']['outname'] + '_Figure3_GC.pdf'
    conf_dict['QCplots']['read_chrom'] = QC_dir + conf_dict['General']['outname'] + '_Figure4_peak_on_chromsome.pdf'
    # conf_dict['QCplots']['peak_over_cis'] = QC_dir + conf_dict['General']['outname'] + '_Distribution_of_Genome_and_ChIP_regions_over_cis_regulatory_element.pdf'
    conf_dict['QCplots']['peak_dis'] = QC_dir + conf_dict['General']['outname'] + '_Figure5_peak_distribution.pdf'
    conf_dict['QCplots']['gene_cover'] = QC_dir + conf_dict['General']['outname'] + '_Figure6_GeneCover.pdf'
    conf_dict['QCplots']['read_dis'] = QC_dir + conf_dict['General']['outname'] + '_Figure7_readsDistribution.pdf'
    conf_dict["results"]["totalPeaks"] = QC_dir + conf_dict['General']['outname'] + "_peaks.bed"

    cmd = "Rscript %s %s"%(conf_dict['rscript']+'DrChIP_readsbulk_QC.r',conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    cmd = "Rscript %s %s %s"%(conf_dict['rscript']+'DrChIP_reads_distribution.r',conf_dict['General']['outputdirectory']+'mapping/'+conf_dict['General']['outname']+'.readsnumber.txt',conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    
    bulkqctime = time.time() - t
    Log("time for quality control: %s"%(bulkqctime),logfile)
    Log("Step2 bulk and individual cell QC DONE",logfile)
    return conf_dict
