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

def CeasSelected(conf_dict,logfile):
    ceas_r_info = open(conf_dict['General']['outputdirectory'] + "QC/" + conf_dict["General"]["outname"]+"_ceas.R")
    tmp_out_info = open(conf_dict['General']['outputdirectory'] + "QC/" + conf_dict["General"]["outname"]+"_ceas_eachplot.R","w")
    figure_marker = 0
    for each in ceas_r_info:
        if each.startswith("mtext"):
            figure_marker = 1
        if each.startswith("pdf"):
            pass
        else:    
            if each.startswith("# Chromosomal Distribution"):
                tmp_out_info.write(each)
                tmp_out_info.write("pdf('%s_Figure4_peak_on_chromsome.pdf',height=12,width=12)\n"%conf_dict["General"]["outname"])
            elif each.startswith("# Promoter,Bipromoter,Downstream, Gene and Regions of interest"):
                tmp_out_info.write("dev.off()\n")
                tmp_out_info.write(each)
                tmp_out_info.write("pdf('%s_peak_fraction_on_promoter_bipromoter_downstream_and_gene.pdf',height=12,width=12)\n"%conf_dict["General"]["outname"])
            elif each.startswith("# Distribution of Genome and ChIP regions over cis-regulatory element"):
                tmp_out_info.write("dev.off()\n")
                tmp_out_info.write(each)
                tmp_out_info.write("pdf('%s_Distribution_of_Genome_and_ChIP_regions_over_cis_regulatory_element.pdf',height=12,width=12)\n"%conf_dict["General"]["outname"])
            elif each.startswith("# ChIP regions over the genome"):
                tmp_out_info.write("dev.off()\n")
                tmp_out_info.write(each)
                tmp_out_info.write("pdf('%s_Figure5_peak_distribution.pdf',height=12,width=12)\n"%conf_dict["General"]["outname"])
            elif each.startswith("par(mar=c(4, 4, 5, 3.8),oma=c(4, 2, 4, 2))") and figure_marker == 1:
                tmp_out_info.write("dev.off()\n")
                tmp_out_info.write("pdf('%s_Figure6_GeneCover.pdf',height=10,width=10)\n"%conf_dict["General"]["outname"])
                tmp_out_info.write(each)
                figure_marker = 0
            elif each.startswith("layout(matrix(c(1, 2, 3, 3, 4, 5), 3, 2, byrow = TRUE),widths=c(1, 1),heights=c(1, 1, 1))"):
                tmp_out_info.write("layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE),widths=c(1, 1),heights=c(1, 1, 1))\n")
            elif each.startswith("par(mfrow=c(3, 2)"):
                tmp_out_info.write("dev.off()\n")
                tmp_out_info.write("pdf('%s_averageProfile_on_exon_intron_region.pdf',height=12,width=12)\n"%conf_dict["General"]["outname"])
                tmp_out_info.write(each)
            else:
                tmp_out_info.write(each)
    ceas_r_info.close()
    tmp_out_info.close()
    os.chdir(conf_dict['General']['outputdirectory'] + "QC/")                
    LogCommand("Rscript "+conf_dict["General"]["outname"]+"_ceas_eachplot.R",logfile)

def Step4Summary(conf_dict,logfile):
    '''

    '''
    # start
    # create section for 
    CeasSelected(conf_dict,logfile)
    
    Log('Step4: summary',logfile)
    Log('copy results',logfile)
    summarydir = conf_dict['General']['outputdirectory'] + 'summary/'
    CreateDirectory(summarydir)
    os.chdir(summarydir)
    
    plot_folder = summarydir + "plots/"
    CreateDirectory(plot_folder)
    os.chdir(plot_folder)
    ### collect results 
    for i in conf_dict['QCplots']:
        if os.path.isfile(conf_dict['QCplots'][i]):
            #realname
            cmd = 'cp %s .'%conf_dict['QCplots'][i]
            LogCommand(cmd,logfile)

    result_folder = summarydir + "results/"
    CreateDirectory(result_folder)
    os.chdir(result_folder)
    for i in conf_dict['results']:
        cmd = 'cp %s .'%conf_dict['results'][i]
        LogCommand(cmd,logfile)

    os.chdir(summarydir)

    Log('generate qc documents',logfile)
    ### initiate 
    QCdoc = """\documentclass[11pt,a4paper]{article}
\usepackage{tabularx}
\usepackage[english]{babel}
\usepackage{array}
\usepackage{graphicx}
\usepackage{color}
\DeclareGraphicsExtensions{.eps,.png,.pdf,.ps}
\\begin{document}
\\title{Dr.seq 2.0 QC and Analysis Summary Report: %s}

\\vspace{-1cm}
\maketitle
\\tableofcontents
\\newpage
\\newpage
\section{Data description}
\\begin{quotation}
Table 1 mainly describe the input file and mapping and analysis parameters.
\end{quotation}
\\begin{table}[h]
\caption{Data description}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }

"""%(LatexFormat(conf_dict['General']['outname']))
    
    ### table1 prepare parameter
    if int(conf_dict['Step1_Mapping']['q30filter']) == 1:
        q30filter = "True"
    else:
        q30filter = "False"
          
    QCdoc += """      
\hline
parameter & value  \\\\
\hline
output name & %s \\\\
\hline
barcode file(file name only) & %s \\\\
\hline
reads file(file name only) & %s \\\\
\hline
reads file format & %s  \\\\
\hline
cell barcode range in read1 &  %s \\\\
\hline
cell barcode range in read2 &  %s \\\\
\hline
mapping software & %s \\\\
\hline
Q30 filter mapped reads & %s \\\\
\hline
maximum fragment length & %s \\\\
\hline
trim bases from 5'/left end of reads & %s \\\\
\hline
threshold for macs14 peak calling & %s \\\\
\hline
size of peak extension & %s \\\\
\hline
the max number of input barcode & %s \\\\
\hline
\end{tabularx}
\end{table}
"""%(LatexFormat(conf_dict['General']['outname']),
     LatexFormat(conf_dict['General']['fastq_1'].split("/")[-1]),
     LatexFormat(conf_dict['General']['fastq_2'].split("/")[-1]),
     conf_dict['General']['format'].upper(),
     conf_dict['General']['barcode_range_1'],
     conf_dict['General']['barcode_range_2'],
     conf_dict['Step1_Mapping']['mapping_software'],
     q30filter,
     conf_dict['Step1_Mapping']['x'],
     conf_dict['Step1_Mapping']['trim5'],
     conf_dict['Step2_QC']['peak_calling_cutoff'],
     conf_dict['Step3_CellClustering']['peak_extend_size'],
     conf_dict['Step3_CellClustering']['max_barcode_num']
     )
    
    ### bulk QC
    QCdoc += """
\\newpage
\\newpage
\section{Reads level QC}
In the reads level QC step we measured the quality of sequencing reads, including nucleotide quality and composition. In the reads level QC step and Bulk-cell level QC step we randomly sampled down total reads to 5 million and used a published package called ``RseQC" for reference.(Wang, L., Wang, S. and Li, W. (2012) )
\\newpage
\\newpage
\subsection{Reads quality}
\\begin{quotation}
Reads quality is one of the basic reads level quality control methods. We plotted the distribution of a widely used Phred Quality Score at every position of sequence to measure the basic sequence quality of your data. Phred Quality Score was calculate by a python function $ord(Q) - 33$. Color in the heatmap represented frequency of this quality score observed at this position. Red represented higher frequency while blue was lower frequency. You may observe a decreasing of quality near the 3'end of sequence because of general degradation of quality over the duration of long runs. If the decreasing of quality influence the mappability (see ``Bulk-cell level QC") then the common remedy is to perform quality trimming where reads are truncated based on their average quality or you can trim serveal base pair near 3'end directly. If it doesn't help, you may consider your Drop-ChIP data poor quality. 
\end{quotation}
\\begin{figure}[h]
        \caption{Reads quality} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

\\newpage
\\newpage
\subsection{Reads nucleotide composition}
\\begin{quotation}
We assess the nucleotide composition bias of a sample. The proportion of four different nucleotides was calculated at each position of reads. Theoretically four nucleotides had similar proportion at each position of reads. You may observe higher A/T count at 3'end of reads because of the 3'end polyA tail generated in sequencing cDNA libaray, otherwise the A/T count should be closer to C/G count. In any case, you should observe a stable pattern at least in the 3'end of reads. Spikes (un-stable pattern) which occur in the middle or tail of the reads indicate low sequence quality. You can trim serveral un-stable bases from the 3'end if low mappability (see ``Bulk-cell level QC") is also observed. If it doesn't help, you may consider your Drop-ChIP data poor quality. Note that t
he A/T vs G/C content can greatly vary from species to species. 
\end{quotation}
\\begin{figure}[h]
        \caption{Reads nucleotide composition} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

\\newpage
\\newpage
\subsection{Reads GC content}
\\begin{quotation}
Distribution of GC content of each read. This module measures the general quality of the library. If the distribution looks different from a single bell (too sharp or too broad) then there may be a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example), which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species. If you observe sharp peak or broder peak and also observe low mappability (see ``Bulk-cell level QC"), you may consider your Drop-ChIP data poor quality.
\end{quotation}
\\begin{figure}[h]
        \caption{Reads GC content} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%((conf_dict['QCplots']['read_qul'].split("/")[-1]),
     (conf_dict['QCplots']['read_nvc'].split("/")[-1]),
     (conf_dict['QCplots']['read_gc'].split("/")[-1])
    )

    QCdoc += """
\\newpage
\\newpage
\section{Bulk-cell level QC}
In the bulk-cell level QC step we measured the performance of total Drop-ChIP reads. In this step we did't separate reads, just like treated the sample as bulk ChIP-seq sample.
\\newpage
\\newpage
\subsection{Reads alignment summary}
\\begin{quotation}
The following table shows reads number after each filter strategy and mapped reads of final selected reads. It measures the general sequencing quality. Low mappability indicates poor sequence quality(see ``Reads level QC") or library quality(caused by contaminant). In summary, if the percentage of ``total mapped reads" is less than 5\\%%, users may consider reconstruct your library(redo the experiment), but first you should make sure you already trim the adapter and map your reads to the corresponded species(genome version). Mappable reads was after Q30 filtering if Q30 filter function was turned on.\\\\
\end{quotation}
\\begin{table}[h]
\caption{Reads alignment summary}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|X| }
    
\hline
genomic region(Category) &  reads number \\\\
\hline
total number of read pairs & %s \\\\
\hline
number of reads pairs after barcode filtered &  %s (%s\\%%)* \\\\
\hline
number of mapped reads & %s (%s\\%%)* \\\\
\hline
number of reads after length filtered & %s (%s\\%%)* \\\\
\hline

\end{tabularx}
\end{table}
"""%(NumberFormat(str(conf_dict['Step2_QC']['total_reads_pair_N'])),
     NumberFormat(str(conf_dict['Step2_QC']['bc_filter_reads_pair_N'])),
     str( round(100*int(conf_dict['Step2_QC']['bc_filter_reads_pair_N'])*1.0/int(conf_dict['Step2_QC']['total_reads_pair_N']), 2)),
     NumberFormat(str(conf_dict['Step2_QC']['mapped_reads_pair_N'])),
     str( round(100*int(conf_dict['Step2_QC']['mapped_reads_pair_N'])*1.0/int(conf_dict['Step2_QC']['total_reads_pair_N']), 2)),
     NumberFormat(str(conf_dict['Step2_QC']['final_reads_pair_N'])),
     str( round(100*int(conf_dict['Step2_QC']['final_reads_pair_N'])*1.0/int(conf_dict['Step2_QC']['total_reads_pair_N']), 2)))
     ### reads on chromsome
    QCdoc += """
\\newpage
\\newpage
\subsection{Chromosomal Distribution of ChIP Regions}
\\begin{quotation}
The blue bars represent the percentages of the whole tiled or mappable regions in the chromosomes (genome background) and the red bars showed the percentages of the whole ChIP. These percentages are also marked right next to the bars. P-values for the significance of the relative enrichment of ChIP regions with respect to the gnome background are shown in parentheses next to the percentages of the red bars.
\end{quotation}
\\begin{figure}[h]
        \caption{Chromosomal Distribution of ChIP Regions} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%((conf_dict['QCplots']['read_chrom'].split("/")[-1]))

### peak distribution
    QCdoc += """
\\newpage
\\newpage
\subsection{Peaks over Chromosomes}
\\begin{quotation}
Barplot show ChIP regions distributed over the genome along with their scores or peak heights. The line graph on the top left corner illustrates the distribution of peak heights (or scores). The red bars in the main plot ChIP regions in the input BED file. The x-axis of the main plot represents the actual chromosome sizes.
\end{quotation}
\\begin{figure}[h]
        \caption{Peaks over Chromosomes} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%((conf_dict['QCplots']['peak_dis'].split("/")[-1]))

### average profile on different genome regions
    QCdoc += """
\\newpage
\\newpage
\subsection{average profile on different genome regions}
\\begin{quotation}
Average profiling within/near important genomic features. The panels on the first row display the average ChIP enrichment signals around TSS and TTS of genes, respectively. The bottom panel represents the average ChIP signals on the meta-gene of 3 kb.
\end{quotation}
\\begin{figure}[h]
        \caption{average profile on different genome regions} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%((conf_dict['QCplots']['gene_cover'].split("/")[-1]))

    QCdoc += """
\\newpage
\\newpage
\section{Individual-cell level QC}
In this step we focused on the quality of individual cell and distinguishing cell reads from background reads
\\newpage
\\newpage
\subsection{Reads distribution}
\\begin{quotation}
Drop-ChIP technology has an innate advantage of detecting individual cell reads and background reads due to the barcode information. This module displays the distribution of reads number in each cell and helps to discard barcodes with high rate of background reads (which usually caused by empty cell barcodes and ambient sequence). We plot the distribution of reads number in each cell barcode (though most of cell barcodes don't contain cells, they still have reads) and observed a bimodal distribution of reads number. The red line show the reads distribution. 
\end{quotation}
\\begin{figure}[h]
        \caption{Reads distribution} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(conf_dict['QCplots']['read_dis'].split("/")[-1])

    QCdoc += """
\\newpage
\\newpage
\section{Cell-clustering level QC}
This step composed by h-clustering based on macs14 peaks.
\\newpage
\\newpage
\subsection{Cell clustering}
\\begin{quotation}
We conducted a h-cluster based on macs14 peaks to measure sample's ability to be separated to different cell subtypes. 
\end{quotation}
\\begin{figure}[h]
        \caption{h-clustering based on peak} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(conf_dict['QCplots']['cell_clustering'].split("/")[-1])
    if os.path.isfile(conf_dict['QCplots']['silhouette']):
        QCdoc += """
\\newpage
\\newpage
\subsection{Silhouette of clustering}
\\begin{quotation}
Silhouette method is used to interprate and validate the consistency within clusters defined in previous steps. A poor Silhouette (e.g. average si $<$ 0.2 ) score indicate that the experiments(if not properly done) may not separate well the subpopulations of cells. If most of your clusters have poor Silhouette score, it may indicate a poor quality of your experiments. 
\end{quotation}
\\begin{figure}[h]
        \caption{Silhouette score for clustered STAMPs} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
 
"""%(conf_dict['QCplots']['silhouette'].split("/")[-1])
    QCdoc += """
\\newpage
\\newpage
\subsection{Clustering heatmap}
\\begin{quotation}
Cell Clustering tree and peak region in each cell. The upper panel represents the hieratical clustering results based on each single cell. The second panel with different colors represents decision of cell clustering. The bottom two panels (heatmap and color bar) represent the "combined peaks" occupancy of each single cell.
\end{quotation}
\\begin{figure}[h]
        \caption{h-clustering heatmap} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(conf_dict["QCplots"]['heatmap'].split("/")[-1])
    QCdoc += """
\\newpage
\\newpage
\subsection{ideogram}
\\begin{quotation}
Cluster specific regions were show in each chromsome.
\end{quotation}
\\begin{figure}[h]
        \caption{ideogram of cluster specific regions} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(conf_dict["QCplots"]['ideogram'].split("/")[-1])


    QCdoc += """
\\newpage
\\newpage
\section{Output list}
\\begin{quotation}
All output files were described in the following table
\end{quotation}
\\begin{table}[h]
\caption{output list}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }
\hline
description & filename \\\\
\hline
peak location matrix for each cell & %s  \\\\
\hline
cells in each cluster & %s  \\\\
\hline
cell type specific peaks per each cluster & %s  \\\\
\hline
cell clustering results with Silhouette Score & %s  \\\\
\hline
combined peaks & %s  \\\\
"""%(LatexFormat(conf_dict["results"]['peak_matrix'].split("/")[-1]),
    LatexFormat(conf_dict["results"]['clusterCells'].split("/")[-1]),
    LatexFormat(conf_dict["results"]['specificPeak'].split("/")[-1]),
    LatexFormat(conf_dict["results"]['cluster_with_silhouette_score'].split("/")[-1]),
    LatexFormat(conf_dict["results"]['totalPeaks'].split("/")[-1]))
    QCdoc += """
\hline
summary QC report & %s \\\\
\hline
\end{tabularx}
\end{table} 
\end{document} 
"""%(LatexFormat(conf_dict['General']['outname'])+"\_summary.pdf")

    os.chdir(plot_folder)

    latexfile = conf_dict['General']['outname'] + '_summary.tex'
    outf = open(latexfile,'w')
    outf.write(QCdoc)
    outf.close()
    cmd = "pdflatex %s"%(latexfile)
    cmd2 = 'cp %s %s'%(conf_dict['General']['outname'] + '_summary.pdf',summarydir)
    if conf_dict['General']['latex'] == 1:
        LogCommand(cmd,logfile)
        LogCommand(cmd,logfile)
        LogCommand(cmd2,logfile)
        for files in os.listdir(plot_folder):
            if os.path.isfile(files) and files[-12:-4] == "_summary":
                if not files[-4:] in ['.tex','.pdf',',png','.txt']:
                    cmd = "rm %s"%(files)
                    LogCommand(cmd,logfile)
        Log('pdflatex was detected in default PATH, generate summary report %s'%('summary/'+conf_dict['General']['outname'] + '_summary.pdf'),logfile)
    else:
        Log('pdflatex was not detected in default PATH, generate summary report .tex file in summary/plots folder, you can move the whole summary/plots/ folder to the environment with pdflatex installed and run cmd in the plots/ folder: "pdflatex %s"'%(conf_dict['General']['outname'] + '_summary.tex'),logfile)
           
    Log('Step5 summary DONE, check %s for final outputs'%(summarydir),logfile)


    return conf_dict









