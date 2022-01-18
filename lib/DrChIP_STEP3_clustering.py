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
import DrSeq2pipe
from DrSeq2pipe.Utility import *

# ================================================================================
# function 
# ================================================================================
def MartrixGenerated(overlap_file,output_name,max_barcode_num):
    overlap_info = open(overlap_file)
    out_info = open(output_name+"_signal.txt","w")
    out_info.write("\tC%s\n"%(str(range(1,int(max_barcode_num+1)))[1:-1].replace(", ","\tC")))
    out_dic = {}
    for each in overlap_info:
        each = each.strip().split()
        region = each[0]
        cell = int(each[1])
        if out_dic.has_key(region):
            out_dic[region][cell-1] += 1
        else:
            out_dic[region] = [0]*int(max_barcode_num)
            out_dic[region][cell-1] += 1
    for each in out_dic:
        out_info.write("%s\t%s\n"%(each,str(out_dic[each])[1:-1].replace(",","\t")))
    overlap_info.close()
    out_info.close()

def Step3Clustering(conf_dict,logfile):
    '''
    analysis part
    Clustering
    '''   
    # start
    # create section for 
    t = time.time()
    Log('Step3: analysis',logfile)
    Log('Cell clustering with own script, based on peaks',logfile)

    analysis_dir = conf_dict['General']['outputdirectory'] + 'analysis/'
    CreateDirectory(analysis_dir)
    os.chdir(analysis_dir)
    
    cmd = """awk -v "OFS=\t" '{print ($1,$2,$3,$1"_"$2"_"$3,$5)}' %sQC/%s_peaks.bed > %s_peaks_rename.bed"""%(conf_dict['General']['outputdirectory'],conf_dict['General']["outname"],conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    cmd = "intersectBed -wo -a %s_peaks_rename.bed -b %smapping/%s.selected.bed | cut -f 4,9 > %s_overlap.txt"%(conf_dict['General']['outname'],conf_dict['General']['outputdirectory'],conf_dict['General']['outname'],conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    # time.sleep(300)
    MartrixGenerated("%s_overlap.txt"%conf_dict['General']['outname'],conf_dict['General']['outname'],int(conf_dict['Step3_CellClustering']['max_barcode_num']))    
    cmd = "intersectBed -wo -a %s_peaks_rename.bed -b %smapping/%s.selected.bed | cut -f 1,2,3,9 | sort | uniq > %s_peak_location.txt"%(conf_dict['General']['outname'],conf_dict['General']['outputdirectory'],conf_dict['General']['outname'],conf_dict['General']['outname'])
    LogCommand(cmd,logfile)
    cmd = "Rscript %s %s %s %s %s %s %s %s"%(conf_dict['rscript']+'DrChIP_cell_clustering.r',conf_dict['General']['outname']+"_signal.txt",conf_dict['General']['outname'],conf_dict['Step3_CellClustering']['cut_height'],conf_dict['Step3_CellClustering']['cell_cutoff'],conf_dict['Step3_CellClustering']['peak_cutoff'],conf_dict['General']['outputdirectory'] + 'mapping/%s.readsnumber.txt'%conf_dict['General']['outname'],conf_dict['Step3_CellClustering']['given_cluster_number'])
    LogCommand(cmd,logfile)
    
    if conf_dict["General"]["species"] == "hs":
        cytoBandfile = conf_dict['rscript'] + "hg19cytoBand.txt.gz"
    elif conf_dict["General"]["species"] == "mm":
        cytoBandfile = conf_dict['rscript'] + "mm9cytoBand.txt.gz"
    else:
        print "Only Homo sapiens(hs) and Mus musculus(mm) are supported for -g option. Please Check your input."
        sys.exit(0)
    peak_per_cell_dic = {}
    peak_per_cell_info = open(conf_dict['General']['outname']+"_peak_location.txt")
    for each in peak_per_cell_info :
        each_list = each.strip().split()
        tkey = "C"+each_list[3]
        if peak_per_cell_dic.has_key(tkey):
            peak_per_cell_dic[tkey] += each
        else:
            peak_per_cell_dic[tkey] = each
    peak_per_cell_info.close()

    k = 0
    specific_peak_files = []
    for each_cluster in os.listdir(analysis_dir):
        if "_cluster" in each_cluster and "_cells.txt" in each_cluster:
            files_info = open(each_cluster)
            tmp_cells = []
            for each_cell in files_info:
                tmp_cells.append(each_cell.strip())
            files_info.close()
            if len(tmp_cells) >= 3:
                k += 1
                files_info = open(each_cluster)
                tmpout_info = open("%s_specific.peak.bed"%each_cluster[:-10],"w")
                for each_cell in files_info:
                    each_cell = each_cell.strip()
                    if peak_per_cell_dic.has_key(each_cell):
                        tmpout_info.write(peak_per_cell_dic[each_cell])
                files_info.close()
                tmpout_info.close()
                specific_peak_files.append("%s_specific.peak.bed"%each_cluster[:-10])
            else:
                Log("%s contains %d cells,discarded."%(each_cluster,len(tmp_cells)),logfile)

    if k > 0:
        if conf_dict["General"]["species"] == "hs":
            cytoBandfile = os.path.join(DrSeq2pipe.__path__[0], "Rscript/hg19cytoBand.txt.gz")
        elif conf_dict["General"]["species"] == "mm":
            cytoBandfile = os.path.join(DrSeq2pipe.__path__[0], "Rscript/mm9cytoBand.txt.gz")
        else:
            print "Only Homo sapiens(hs) and Mus musculus(mm) are supported for -g option. Please Check your input."
            sys.exit(0)
        cmd = "Rscript %s %s %s %s %s %s"%(conf_dict['rscript']+"DrChIP_ideogram_draw.r",cytoBandfile,conf_dict['General']['outname'],str(k),",".join(specific_peak_files),conf_dict["General"]["species"])
        LogCommand(cmd,logfile)
    else:
        LogError("Cluster cell files can't be detected in the analysis directory, perhaps too little informative reads can be mapped.Please Check your input.",logfile)

    Log("Step3 cell clustering DONE",logfile)
    if not conf_dict.has_key("results"):
        conf_dict["results"] = {}
    conf_dict['results']['peak_matrix'] = analysis_dir + conf_dict['General']['outname'] + '_signal.txt'
    conf_dict["results"]["clusterCells"] = analysis_dir + conf_dict['General']['outname'] + "_cluster*_cells.txt"
    conf_dict["results"]["specificPeak"] = analysis_dir + conf_dict['General']['outname'] + "*_specific.peak.bed"
    conf_dict["results"]["cluster_with_silhouette_score"] = analysis_dir + conf_dict['General']['outname'] + "_cluster_with_silhouette_score.txt"
    conf_dict['QCplots']['cell_clustering'] = analysis_dir + conf_dict['General']['outname'] + '_Figure8_cell_clusting.pdf'
    conf_dict['QCplots']['silhouette'] = analysis_dir + conf_dict['General']['outname'] + '_Figure9_silhouetteScore.pdf'
    conf_dict['QCplots']['heatmap'] = analysis_dir + conf_dict['General']['outname'] + '_Figure10_heatmap.png'
    conf_dict['QCplots']['ideogram'] = analysis_dir + conf_dict['General']['outname'] + '_Figure11_ideogram.pdf'
    return conf_dict


