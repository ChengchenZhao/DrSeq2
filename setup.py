#!/usr/bin/env python
"""
Description
Setup script for DrSeq2  -- QC and analysis pipeline for parallel single cell transcriptome data and epigenome data
Copyright (c) 2016 <Chengchen Zhao> <1310780@tongji.edu.cn>
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""
import os
import sys
import subprocess
import platform
from distutils.core import setup, Extension

def setup_ceas():
    setup(name="CEAS-Package",
      version="1.0.2",
      description="CEAS -- Cis-regulatory Element Annotation System Package",
      author='H. Gene Shin',
      author_email='shin@jimmy.harvard.edu',
      url='http://liulab.dfci.harvard.edu/CEAS/',
      package_dir={'CEAS' : 'ceas_lib'},
      packages=['CEAS'],
      scripts=['bin/ceas', 'bin/sitepro', 'bin/gca', 'bin/build_genomeBG'],

      classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Artistic License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Database',
        ],
      )

def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
def compile_bedtools():
    curdir = os.getcwd()
    os.chdir('refpackage/bedtools')
    sp('make 1>/dev/null 2>&1 ')
    sp('chmod 755 *')
    print("bedtools has been installed.")
    os.chdir(curdir)
def check_bedtools():
    checkhandle_1 = sp('which bedtools')
    checkhandle_2 = sp('which bamToBed')
    if checkhandle_1[0].strip() == "" or checkhandle_2[0].strip() == "":
        return 0
    else:
        return 1
# def check_bedGraphToBigWig():
#     checkhandle = sp('which bedGraphToBigWig')
#     if checkhandle[0].strip() == "":
#         return 0
#     else:
#         return 1
def check_samtools():
    checkhandle = sp('which samtools')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1
# def check_bedClip():
#     checkhandle = sp('which bedClip')
#     if checkhandle[0].strip() == "":
#         return 0
#     else:
#         return 1
def check_R():
    checkhandle = sp('which Rscript')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1    
def compile_SIMLR():
    curdir = os.getcwd()
    os.chdir('lib/Rscript/')
    if os.path.isfile("lib/Rscript/projsplx_R.so"):
        sp("rm projsplx_R.so")
    sp("command R CMD SHLIB -c projsplx_R.c")
    os.chdir(curdir)

def main(): 
    if sys.version_info[0] != 2 or sys.version_info[1] < 7:
        print >> sys.stderr, "ERROR: DrSeq2 requires Python 2.7"
        sys.exit()
    has_R = check_R()
    if has_R == 0:
        print >> sys.stderr, "ERROR: DrSeq2 requires R & Rscript under default PATH"
        sys.exit()        
    has_bedtools = check_bedtools()
    try:
        import CEAS.inout as inout
    except ImportError:
        print "CEAS is not detected under default PATH, install CEAS first."
        setup_ceas()
    # R packages
    # print 'Intalling R packages for DrSeq2.'
    # os.system("Rscript install_R_packages.r")
    # print "Necessary R package has been installed."
    
    # different version of bedGraphtoBigwig/bedClip/samtools for different system
    if platform.system() == "Linux":
        script_dir = "bin/linux.x86_64/"
    elif platform.system() == "Darwin":
        script_dir = "bin/macOSX.i386/"
    else:
        script_dir = ""
    
    scripts_list=['bin/DrSeq','bin/DrChIP','bin/ceas','bin/sitepro','bin/gca','bin/build_genomeBG','bin/10x2Dr','bin/MARS2Dr','bin/comCluster','bin/GeMa','bin/ATAC']
    # bedtools
    if has_bedtools == 0:
        print 'bedtools is not detected under default PATH, install bedtools for DrSeq2, may take serval minutes...'
        compile_bedtools()
        scripts_list.append('refpackage/bedtools/bin/bedtools')
        scripts_list.append('refpackage/bedtools/bin/bamToBed')
        scripts_list.append('refpackage/bedtools/bin/intersectBed')

    # if not check_bedGraphToBigWig():
    #     if script_dir == "":
    #         print 'bedGraphToBigWig is required. For convenience we provide the executable file for macOSX and Linux system. But your system is %s'%platform.system()
    #         sys.exit()
    #     else:
    #         scripts_list.append(script_dir+"bedGraphToBigWig")
    if not check_samtools():
        if script_dir == "":
            print 'samtools is required. For convenience we provide the executable file for macOSX and Linux system. But your system is %s. Please ensure samtools installed and try DrSeq2 installation again.'%platform.system()
            sys.exit()
        else:
            scripts_list.append(script_dir+"samtools")
    # if not check_bedClip():
    #     if script_dir == "":
    #         print 'bedClip is required. For convenience we provide the executable file for macOSX and Linux system. But your system is %s'%platform.system()
    #         sys.exit()
    #     else:
    #         scripts_list.append(script_dir+"bedClip")

    setup(name="DrSeq2pipe",
          version="2.2.0",
          description="DrSeq2 : a quality control and analysis pipeline for parallel single cell transcriptome and epigenome data",
          author='Chengchen Zhao',
          author_email='1310780@tongji.edu.cn',
          url='https://github.com/ChengchenZhao/DrSeq2.0.git',
          package_dir={'DrSeq2pipe' : 'lib'},
          packages=['DrSeq2pipe'],
          package_data={'DrSeq2pipe': ['Config/*','Rscript/*','Data/*',]},
          scripts=scripts_list,
          classifiers=[
          'Development Status :: version1.0 finish',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: GPLv3',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          'Topic :: pipeline',
          ],
          requires=[],
    )
        
if __name__ == "__main__":
    print 'Intalling DrSeq2, may take serval minutes...'
    main()
    print 'Installation of DrSeq2 is DONE!'

