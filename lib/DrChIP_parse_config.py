#!/usr/bin/env python
# ------------------------------------
"""
Function declare:

def NewConf(conf_name)
def ReadConf(conf_file)
def MakeConf(fastq_1,fastq_2,barcode_file,barcode_file_range_1,barcode_file_range_2,barcode_range_1,barcode_range_2,outname,geneanno,maptool,mapindex,P,X,trim5,fover,Clean,species,cell_cutoff,peak_cutoff):

"""
# -----------------------------------


import os,sys
import ConfigParser

import DrSeq2pipe

### config template
CONFIG_TEMPLATE = os.path.join(DrSeq2pipe.__path__[0], "Config/DrChIP_template.conf")

### generate a config
def NewConf(conf_name):
    '''
    Generate a config file
    ''' 
    inf = open(CONFIG_TEMPLATE)
    if not conf_name.endswith('.conf'):
        conf_name += '.conf'
    outf = open(conf_name,'w')
    for line in inf:
        outf.write(line)
    outf.close()
    inf.close()

### read config
def ReadConf(conf_file):
    '''
    Read config file and return a dict containing all infomation
    '''
    conf_dict = {}
    cf = ConfigParser.SafeConfigParser()
    cf.read(conf_file)
    for st in cf.sections():
        conf_dict[st]={}
        for item in cf.items(st):
            conf_dict[st][item[0]]=item[1]
    return conf_dict

### generate a config file in simple mode
def MakeConf(fastq_1,fastq_2,barcode_file,barcode_file_range_1,barcode_file_range_2,barcode_range_1,barcode_range_2,outname,geneanno,maptool,mapindex,P,X,trim5,fover,Clean,species,cell_cutoff,peak_cutoff):
    inf = open(CONFIG_TEMPLATE)
    if os.path.isfile(outname+'.conf') and not fover :
        print 'config file "%s" using same name exist , choose other name or add -f to overwrite'%(outname+".conf")
        sys.exit(1)
    outf = open(outname+".conf",'w')
    for line in inf:
        if line.startswith('fastq_1 ='):
            newline = 'fastq_1 = ' + fastq_1 + '\n'
        elif line.startswith('fastq_2 ='):
            newline = 'fastq_2 = ' + fastq_2 + '\n'
        elif line.startswith('barcode_file ='):
            newline = 'barcode_file = ' + str(barcode_file) + '\n'
        elif line.startswith('barcode_file_range_1 ='):
            newline = 'barcode_file_range_1 = ' + str(barcode_file_range_1) + '\n'
        elif line.startswith('barcode_file_range_2 ='):
            newline = 'barcode_file_range_2 = ' + str(barcode_file_range_2) + '\n'
        elif line.startswith('barcode_range_1 ='):
            newline = 'barcode_range_1 = ' + str(barcode_range_1) + '\n'
        elif line.startswith('barcode_range_2 ='):
            newline = 'barcode_range_2 = ' + str(barcode_range_2) + '\n'
        elif line.startswith('outputdirectory ='):
            newline = 'outputdirectory = ' + outname + '\n'
        elif line.startswith('outname ='):
            newline = 'outname = ' + outname + '\n'
        elif line.startswith('mapping_software ='):
            newline = 'mapping_software = ' + maptool + '\n'
        elif line.startswith('gene_annotation ='):
            if geneanno:
                newline = 'gene_annotation = ' + geneanno + '\n'
            else:
                newline = line
        elif line.startswith('mapindex ='):
            if mapindex:
                newline = 'mapindex = ' + mapindex + '\n'
            else:
                newline = line
        elif line.startswith('p ='):
            newline = 'p = ' + str(P) + '\n'
        elif line.startswith('X ='):
            newline = 'X = ' + str(X) + '\n'
        elif line.startswith('trim5 ='):
            newline = 'trim5 = ' + str(trim5) + '\n'
        elif line.startswith('species ='):
            newline = 'species = ' + str(species) + '\n'
        elif line.startswith('cell_cutoff ='):
            newline = 'cell_cutoff = ' + str(cell_cutoff) + '\n'
        elif line.startswith('peak_cutoff ='):
            newline = 'peak_cutoff = ' + str(peak_cutoff) + '\n'
        else:
            newline = line
        outf.write(newline)
    inf.close()
    outf.close()
    return  
