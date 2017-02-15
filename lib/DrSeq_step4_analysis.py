#!/usr/bin/env python

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
                                   CreateDirectory)
# --------------------------
# main 
# --------------------------
def step4_analysis(conf_dict,logfile):
    '''
    analysis part
    mainly Rscript
    dimentional reduction + clustering
    '''   
    # start
    # create section for 
    t = time.time()
    Log('Step4: analysis',logfile)
    Log('dimentional reduction + clustering with own script, based on selected STAMP barcodes',logfile)
    analysisdir = conf_dict['General']['outputdirectory'] + 'analysis/'
    CreateDirectory(analysisdir)
    os.chdir(analysisdir)

    conf_dict['Step4_Analysis']['clusterresult'] = analysisdir + conf_dict['General']['outname']+'_cluster.txt'
    conf_dict['QCplots']['gapstat'] = analysisdir + conf_dict['General']['outname']+'_Figure10_GapStat.pdf'
    conf_dict['QCplots']['cluster'] = analysisdir + conf_dict['General']['outname']+'_Figure11_cluster.pdf'
    conf_dict['QCplots']['silhouette'] = analysisdir + conf_dict['General']['outname']+'_Figure12_silhouetteScore.pdf'
    conf_dict['QCplots']['umicolor'] = analysisdir + conf_dict['General']['outname']+'_Figure13_totalUMIcolored.pdf'
    conf_dict['QCplots']['itrcolor'] = analysisdir + conf_dict['General']['outname']+'_Figure14_intronRate_colored.pdf'
    conf_dict['results']['cortable'] = analysisdir + conf_dict['General']['outname']+'_correlation_table.txt' 
    conf_dict['results']['features'] = analysisdir + conf_dict['General']['outname']+'_pctablefeatures_clustercell.txt'
        
    if int(conf_dict["Step4_Analysis"]["dimensionreduction_method"]) == 1:
        conf_dict['results']['pctable'] = analysisdir + conf_dict['General']['outname']+'_pctable.txt'
        Log("Using t-SNE to ceducted the dimension.",logfile)
        cmd = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s"%('Rscript',conf_dict['rscript']+'DrSeq_analysis.r',conf_dict['results']['expmatcc'],conf_dict['General']['outname'],conf_dict['Step4_Analysis']['highvarz'],conf_dict['Step4_Analysis']['selectpccumvar'],conf_dict['Step4_Analysis']['rdnumber'],conf_dict['Step4_Analysis']['maxknum'],conf_dict['Step4_Analysis']['pctable'],conf_dict['Step4_Analysis']['cortable'],conf_dict['Step4_Analysis']['clustering_method'],conf_dict['Step4_Analysis']['custom_k'],conf_dict['Step4_Analysis']['custom_d'],conf_dict['Step4_Analysis']['seed'])
        LogCommand(cmd,logfile)
    elif int(conf_dict["Step4_Analysis"]["dimensionreduction_method"]) == 2:
        Log("Using SIMLR to ceducted the dimension.",logfile)
        conf_dict['Step4_Analysis']['pctable'] = 0
        cmd = "%s %s %s %s %s %s %s %s %s %s %s %s"%('Rscript',conf_dict['rscript']+'DrSeq_SIMLR.r',conf_dict['results']['expmatcc'],conf_dict['General']['outname'],conf_dict['Step4_Analysis']['highvarz'],conf_dict['Step4_Analysis']['rdnumber'],conf_dict['Step4_Analysis']['maxknum'],conf_dict['Step4_Analysis']['cortable'],conf_dict['Step4_Analysis']['clustering_method'],conf_dict['Step4_Analysis']['custom_k'],conf_dict['Step4_Analysis']['custom_d'],conf_dict['rscript'])
        LogCommand(cmd,logfile)
    elif int(conf_dict["Step4_Analysis"]["dimensionreduction_method"]) == 3:
        Log("Using PCA to ceducted the dimension.",logfile)
        cmd = "%s %s %s %s %s %s %s %s %s %s %s %s %s"%('Rscript',conf_dict['rscript']+'DrSeq_PCA.r',conf_dict['results']['expmatcc'],conf_dict['General']['outname'],conf_dict['Step4_Analysis']['highvarz'],conf_dict['Step4_Analysis']['rdnumber'],conf_dict['Step4_Analysis']['maxknum'],conf_dict['Step4_Analysis']['pctable'],conf_dict['Step4_Analysis']['cortable'],conf_dict['Step4_Analysis']['clustering_method'],conf_dict['Step4_Analysis']['custom_k'],conf_dict['Step4_Analysis']['custom_d'],conf_dict['rscript'])
        LogCommand(cmd,logfile)
        conf_dict['results']['pctable'] = analysisdir + conf_dict['General']['outname']+'_pctable.txt'
    else:
        Log("You can only choose the method of Dimentional reduction from t-SNE,SIMLR and PCA by now.!",logfile)
    cmd = '%s %s %s %s %s'%('Rscript',conf_dict['rscript']+'DrSeq_post_analysis.r',conf_dict['Step4_Analysis']['clusterresult'],conf_dict['Step2_ExpMat']['qcmatcc'],conf_dict['General']['outname'])


    analysisqctime = time.time() - t
    Log("time for analysis qc: %s"%(analysisqctime),logfile)
    Log("Step4 analysis QC DONE",logfile)
    
    return conf_dict


