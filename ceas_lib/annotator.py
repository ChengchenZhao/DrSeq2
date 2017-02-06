
"""Module Description

Copyright (c) 2008 H. Gene Shin <shin@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  H. Gene Shin
@contact: shin@jimmy.harvard.edu
"""

# ------------------------------------
# Python modules
# ------------------------------------
import sys,time,re,operator,copy,sqlite3,warnings
import itertools
from array import *
import bisect

# ------------------------------------
# My own Python modules
# ------------------------------------
from CEAS.inout import *
from CEAS.tables import *
from CEAS.sampler import *
from CEAS.corelib import *
import CEAS.Prob as Prob

#-------------------------------------
# classes
#-------------------------------------  
class Annotator:
    """Annotator Class
    
    This class annotates a list of genome coordinates and gives a summary of annotation.
    
    1. Annotater.annotate() annotates a list of genome coordinates (locations) based on a gene annotation table (e.g., refSeq provided by UCSC).
    2. Annotatoer.summarize() summarizes the annotation in a table (See tables.py)
    
    """
    
    def __init__(self): 
        """Constructor"""        
        pass
        
    def annotate(self,genome_coordinates=None,gene_table=None,roi=None,prom=(1000, 2000, 3000),biprom=(2500, 5000), down=(1000, 2000, 3000), gene_div=(3,5),quantize=True):
        """Annotate given coordinates based on the given gene table."""
        
        # get the chromsomes of the gene table and genome coordinates
        try:
            chroms_gc=genome_coordinates.keys()
            chroms_gt=gene_table.get_chroms()
            chroms_gt,chroms_gc=set(chroms_gt),set(chroms_gc)
            chroms=chroms_gt.intersection(chroms_gc)
            chroms=list(chroms)
            chroms.sort()
            num_coordinates={}
            num_genes={}
            for chrom in chroms:
                num_coordinates[chrom]=len(genome_coordinates[chrom])
                num_genes[chrom]=len(gene_table[chrom][gene_table[chrom].keys()[0]])
        except AttributeError:
            raise Exception('Genome coordinates and gene table must be given for genome annotation')
        
        
        #initialize with an empty dictionary
        table=AnnotTable()
        #iterate over the chromosomes
        for chrom in chroms: 
            genes=gene_table[chrom]
            num_genes_this_chr=num_genes[chrom]
            coordinates=genome_coordinates[chrom]
            num_coordinates_this_chr=num_coordinates[chrom]
            table.init_table(chrom)
            
            #initialize the promoter distances. This will be used in obtaining bidirectional promoter too.
            prom_dists=[[0,0] for i in xrange(num_coordinates_this_chr)]
        
            # get the nearest genes to set the searching range for promoter and downstream
            nearest_genes=self.find_nearest_genes(genes)
            
            # point begin
            pointerBeg=0
            maxprom = max(prom[-1], biprom[-1])
            maxdown = down[-1]
            for i in xrange(0,num_genes_this_chr):

                # get the strand of the gene
                try:
                    strand=genes['strand'][i]
                except KeyError:
                    raise Exception("'strand' must be included in the gene annotation table for running CEAS")
                
                # get the beginning and end point of search
                # the beginning and end points are the end of the previous gene and the beginning of the next gene.
                beg,end=0,0
                try:
                    if strand=='+':
                        beg=max(genes['txStart'][i]-maxprom, nearest_genes['before'][i])
                        end=min(genes['txEnd'][i]+maxdown, nearest_genes['after'][i])
                    else:
                        beg=max(genes['txStart'][i]-maxdown, nearest_genes['before'][i])
                        end=min(genes['txEnd'][i]+maxprom, nearest_genes['after'][i])
                except KeyError:    # check the gene annotation table has necessary columns
                    raise Exception("'txStart' and 'txEnd' must be included in the gene annotation table for running CEAS")
                
                ### test block-out ###
                # set search index j to the begining point of the last gene. This makes sure that we include isoforms
                j=pointerBeg
                ###

                ### test block-in ###
                #j = bisect.bisect_left(coordinates, beg)
                ###
                if coordinates[j]>end: continue
                
                ### test block-out ###
                # two while loops to detect the annotation start coordinate for the current gene.
                while j>0 and coordinates[j]>=beg: 
                    j-=1  
                while j<num_coordinates_this_chr and coordinates[j]<beg:
                    if j>=table.size(chrom)[0]:
                        table.add_row(chrom,[coordinates[j]]+[0]*table.get_column_num())
                    j+=1
                ###

                # if get to the end of chromosome, then break
                if j==num_coordinates_this_chr: break
                
                ### test block-out
                # save the current start point for the next gene
                pointerBeg=j
                ###

                # otherwise, get the annotations of the probes related with the current gene
                while j<num_coordinates_this_chr and (coordinates[j]>=beg and coordinates[j]<=end):
                    # get the annotation and update the entire annotation table
                    single_annot=self.annotate_single(coordinates[j],strand,genes['txStart'][i],genes['txEnd'][i],\
                                                       genes['cdsStart'][i],genes['cdsEnd'][i],genes['exonStarts'][i],genes['exonEnds'][i],prom,down,gene_div,biprom)
                    self.update_annot_table(table,single_annot,coordinates,prom_dists,chrom,j,biprom)
                    j+=1
            
            # quantize promoter, bipromoter and downstream
            if quantize:
                table[chrom]['promoter']=ordinate2(table[chrom]['promoter'],prom)
                table[chrom]['bipromoter']=ordinate2(table[chrom]['bipromoter'],biprom)
                table[chrom]['downstream']=ordinate2(table[chrom]['downstream'],down)
        
        if roi:
            roichroms = roi.get_chroms()
            for chrom in chroms:
                table[chrom]['roi']=[0]*len(genome_coordinates[chrom])
                if chrom in roichroms:
                    self.do_annotation_roi(genome_coordinates[chrom], table[chrom], roi[chrom])
       
        return table
    
    
                  
    def annotate_single(self,coordinate,strand,txStart,txEnd,cdsStart,cdsEnd,exonStarts,exonEnds,prom,down,gene_div,biprom):
        """Annotate a single genome coordinate
        
        Parameters:
        1. coordinate: a list of genome locations to annotate
        2. strand: the strand (+/-) of the gene
        3. txStart: transcription start
        4. txEnd: transcription end
        5. cdsStart: translation start
        6. cdsEnd: translation end
        7. exonStarts: exon start locations
        8. exonEnds: exon end locations
        9. prom: promoter lengths (e.g., (1000, 2000, 3000))
        10. down: downstream lengths (e.g., (1000, 2000, 3000))
        11. gene: the number of divisions of a gene (eg, (3,5))
        
        """
                               
        # container of the annotation for a single location. 
        single_annot=None
        
        maxprom = max(prom[-1], biprom[-1])
        maxdown = down[-1]
        # get the annotation for the location
        # + strand
        if strand=='+':
            # promoter
            if coordinate<txStart and coordinate>=txStart-maxprom: # txStart-promo <= prob < txStart
                single_annot=['promoter',txStart-coordinate]
            # downstream     
            elif coordinate>=txEnd and coordinate<txEnd+maxdown:# txEnd <= prob < txEnd+down
                single_annot=['downstream',coordinate-txEnd]
            # gene
            elif coordinate>=txStart and coordinate<txEnd:    
                isExon=self.exonify(coordinate,exonStarts,exonEnds) 
                # exon
                if isExon != None:
                    # coding exon
                    if coordinate>=cdsStart and coordinate<cdsEnd: single_annot=['gene',3,self.get_rel_loc_cds(coordinate,strand,exonStarts,exonEnds,isExon,gene_div)]
                    # 5'UTR
                    elif coordinate<cdsStart: single_annot=['gene',1]
                    # 3'UTR
                    else: single_annot=['gene',2]
                # intron
                else: 
                    single_annot=['gene',4]
                # relative location within the gene
                single_annot.insert(2,self.get_rel_loc(coordinate,strand,txStart,txEnd,gene_div))
        # - strand
        else: 
            # promoter
            if coordinate>=txEnd and coordinate<txEnd+maxprom:  
                single_annot=['promoter',txEnd-coordinate]
            # downstream
            elif coordinate<txStart and coordinate>=txStart-maxdown:
                single_annot=['downstream',txStart-coordinate]
            # gene
            elif coordinate>=txStart and coordinate<txEnd:
                isExon=self.exonify(coordinate,exonStarts,exonEnds) #then in an exon
                # exon
                if isExon != None:
                    # coding exon
                    if coordinate>=cdsStart and coordinate<cdsEnd: single_annot=['gene',3,self.get_rel_loc_cds(coordinate,strand,exonStarts,exonEnds,isExon,gene_div)]
                    # 5'UTR
                    elif coordinate>=cdsEnd: single_annot=['gene',1]
                    # 3'UTR
                    else: single_annot=['gene',2]
                # intron
                else: 
                    single_annot=['gene',4]
                # relative location within the gene
                single_annot.insert(2,self.get_rel_loc(coordinate,strand,txStart,txEnd,gene_div))
        
        return single_annot
    
              
    def exonify(self,coordinate,exonStarts,exonEnds):
        """Return which exon a coordinate falls in"""
            
        hit=None
        for i in xrange(0, len(exonStarts)):
            try:
                if coordinate>=exonStarts[i] and coordinate<exonEnds[i]: hit=i
            except IndexError:
                raise Exception("'exonStarts' and 'exonEnds' must have the same number of elments")
        return hit
    
    def get_rel_loc(self,coordinate,strand,txStart,txEnd,gene_div):
        """Return the relative location within the gene"""

        lengths=[(i,int(round(1.0*(txEnd-txStart+1)/i))) for i in gene_div]
        if strand=='+':
            return [min((coordinate-txStart+1)/length+1,i) for i,length in lengths]
        else:
            return [min((txEnd-coordinate+1)/length+1,i) for i,length in lengths]
        
    
    def get_rel_loc_cds(self,coordinate,strand,exonStarts,exonEnds,isExon,gene_div):
        """Return the relative locaiton within the CDS"""

        array_subtractor=lambda array1,array2: array('l',[array1[i]-array2[i] for i in range(0,len(array1))])
        cdsLength=sum(array_subtractor(exonEnds,exonStarts))
        lengths=[(i,int(round(1.0*cdsLength/i))) for i in gene_div]
        
        if strand=='+':
            return [min((coordinate-exonStarts[isExon]+sum(array_subtractor(exonEnds[:isExon],exonStarts[:isExon])))/length+1,i) for i,length in lengths]
        else:
            return [min((exonEnds[isExon]-coordinate+sum(array_subtractor(exonEnds[isExon+1:],exonStarts[isExon+1:])))/length+1,i) for i,length in lengths]        
    
    
    def update_annot_table(self,table,single_annot,coordinates,prom_dists,chrom,index_coord,biprom):
        """Update the annotation table by comparing the current annotation (single_annot) and the existing annotation in the table"""
        
        categories=['promoter','bipromoter','downstream','gene']
        # If the current annotation is NoneType, which means that the curren genome coordinate does not belong to any of
        # promoters,downstreams,or genes, just add an empty annotation to the end of the annotation table
   
        maxbiprom = biprom[-1]
        try:
            category=single_annot[0]
        except TypeError:        # when single_annot=None
            if index_coord>=table.size(chrom)[0]:
                table.add_row(chrom,[coordinates[index_coord]]+[0]*table.get_column_num())
            return
                   
        # When an annotation does not exist at the current genome coordinate, append the new annotation
        # and update
        try:
            original=table[chrom][category][index_coord]
        except IndexError:
                table.add_row(chrom,[coordinates[index_coord]]+[0]*table.get_column_num())
                original=table[chrom][category][index_coord]
    
        # +/- promoter distances
        prom_dist=prom_dists[index_coord]
        
        # when an annotation already exists at the current genome coordinate
#        scaler=lambda dist,maxrange,scale:int(1.0*max(1,dist-1)*scale/maxrange)+1
        new=single_annot[1]
        
        # promoter and bipromoter
        if category=='promoter':
    
            # promoter
            # update + and - distances
            if new < 0:
                if prom_dist[0]==0: 
                    prom_dist[0]=new
                else: 
                    prom_dist[0]=max(new,prom_dist[0])
            elif new > 0:
                if prom_dist[1]==0: 
                    prom_dist[1]=new
                else: 
                    prom_dist[1]=min(new,prom_dist[1])
                        
            # update the promoter annnotation of the current probe
            if original==0:
                if prom_dist[0]==0 or prom_dist[1]==0:
                    fin_prom=abs(sum(prom_dist))
                else: 
                    fin_prom=min(abs(prom_dist[0]),prom_dist[1])
            else: 
                if prom_dist[0]==0 and prom_dist[1]==0:
                    fin_prom=original
                elif prom_dist[0]==0 and prom_dist[1]!=0:
                    fin_prom=min(prom_dist[1],original)
                elif prom_dist[0]!=0 and prom_dist[1]==0:
                    fin_prom=min(abs(prom_dist[0]),original)
                else:
                    fin_prom=min([abs(prom_dist[0]),prom_dist[1],original])
            #update the table with the shortest distance from TSS
            table[chrom][category][index_coord]=fin_prom
            
            # bidirectional promoter
            # if + and - distances have non-zero values, that is bidirectional promoter
            if prom_dist[0]*prom_dist[1] != 0:    
                new_biprom=abs(prom_dist[0])+prom_dist[1]

                #if new_biprom<=cutoff:
                if new_biprom<=maxbiprom:
                    # original bipromoter value
                    ori_biprom=table[chrom]['bipromoter'][index_coord]
                    # update the bipromoter
                    if ori_biprom==0: fin_biprom=new_biprom
                    else: fin_biprom=min(ori_biprom,new_biprom)
                    table[chrom]['bipromoter'][index_coord]=fin_biprom
        elif category=='downstream':
            if original==0: fin_down=new
            else: fin_down=min(new,original)
            table[chrom][category][index_coord]=fin_down
        elif category=='gene':
            if original==0: fin_gene=new
            else: fin_gene=min(new,original)
            table[chrom][category][index_coord]=fin_gene
            table[chrom]['rel_loc'][index_coord]=copy.deepcopy(single_annot[2])
            try:
                table[chrom]['rel_loc_cds'][index_coord]=copy.deepcopy(single_annot[3])
            except IndexError:
                pass
    
    
    def find_nearest_genes(self,genes):
        """Given a list of genes, find their nearest genes ahead of them
        
        Parameters:
        1. genes: a dictionary of gene annotation table. If GeneTable object (see inout.py) is used, simply GeneTable[chrom].
        
        Return:
        nearest_genes: {'before':[n1,n2,...],'after':[m1,m2,...]}. The ith elements under 'txStart' and 'txEnd' represent 
        the end and start of the nearest genes before and after the ith gene.
        
        """

        nearest_genes={'before':[],'after':[]}
        num_genes=len(genes[genes.keys()[0]])
        for i in xrange(num_genes):
            j=i+1
            while j<num_genes and genes['txEnd'][i]>=genes['txStart'][j]:
                j+=1
                
            j=min(num_genes-1,j)
            if genes['txEnd'][i]>=genes['txStart'][j]: 
                nearest_genes['after'].append(genes['txEnd'][num_genes-1]+1000000)
            else:
                nearest_genes['after'].append(genes['txStart'][j])
            
            # find the nearest gene before the current gene
            j=i-1
            while j>=0 and genes['txStart'][i]<=genes['txEnd'][j]:
                j-=1
                
            j=max(0,j)
            if genes['txStart'][i]<=genes['txEnd'][j]:
                nearest_genes['before'].append(0)
            else:
                nearest_genes['before'].append(genes['txEnd'][j])
        
        return nearest_genes


    def do_annotation_roi(self,coordinates,table_of_chrom,roi_of_chrom):
        """Perform annotation for the non-coding regions"""
        
        init = 0
        union = union_intervals([roi_of_chrom['start'],roi_of_chrom['end']])
        for ncstart,ncend in itertools.izip(union[0], union[1]):
        #for ncstart,ncend in itertools.izip(roi_of_chrom['start'],roi_of_chrom['end']):
            start,end=where(ncstart,ncend,coordinates[init:])
            table_of_chrom['roi'][init+start:init+end]=[1]*(end-start)
            init+=end
            
    def summarize(self,table):
        """Provides a summary of the annotation.
        
        Parameters:
        table: an AnnotTable object produced by Annotater.annotate()
        Output:
        summary: a Summary object (see tables.py), which contains the summary of the annotation
        """
        
        # if empty annotation table, just return without doing any.
        # empty summary table
        summary=Summary()
        p=P()
        try:
            chroms=table.get_chroms()
            if not chroms: raise ValueError
        except (AttributeError,ValueError):
            return summary
                    
        # obtain a summary statistics
        summary.init_table('whole')
        for chrom in chroms:
            # get the summary of a single chromosome
            self.get_summary(table,summary,chrom)
            # integrate it into the 'whole'
            self.integrate_summaries(summary,chrom)
            # get p(probablities)
        
        p.init_table('whole')
        for chrom in chroms:
            self.get_p(summary,p,chrom)
        # get the probabilities of the whole genome
        self.get_p(summary,p,'whole')
        
        return summary,p
            
    def get_summary(self,table,summary,chrom):   
        """Get the summary of a single chromosome"""
        
        length=table.size(chrom)[0]
        summary.init_table(chrom)
        summary[chrom]['Ns']=length
                    
        for i in xrange(0,length):
            # get some sort of histograms
            
            if table[chrom]['promoter'][i]!=0: summary[chrom]['promoter'][abs(table[chrom]['promoter'][i])-1]+=1
            if table[chrom]['bipromoter'][i]!=0: summary[chrom]['bipromoter'][table[chrom]['bipromoter'][i]-1]+=1
            if table[chrom]['downstream'][i]!=0: summary[chrom]['downstream'][table[chrom]['downstream'][i]-1]+=1
            
            if table[chrom]['gene'][i]!=0: summary[chrom]['gene'][table[chrom]['gene'][i]-1]+=1
            
            if table[chrom]['rel_loc'][i]!=[0,0]: 
                summary[chrom]['rel_loc'][0][table[chrom]['rel_loc'][i][0]-1]+=1
                summary[chrom]['rel_loc'][1][table[chrom]['rel_loc'][i][1]-1]+=1
            if table[chrom]['rel_loc_cds'][i]!=[0,0]:
                summary[chrom]['rel_loc_cds'][0][table[chrom]['rel_loc_cds'][i][0]-1]+=1
                summary[chrom]['rel_loc_cds'][1][table[chrom]['rel_loc_cds'][i][1]-1]+=1
            try:
                if table[chrom]['roi'][i]!=0: summary[chrom]['roi']+=1
            except KeyError:
                pass
        
        # the last gene element is all
        summary[chrom]['gene'][-1]=sum(summary[chrom]['gene'][:-1]) 
        # get the cumulative sums
        summary[chrom]['promoter']=[sum(summary[chrom]['promoter'][:i]) for i in range(1,len(summary[chrom]['promoter'])+1)]
        summary[chrom]['bipromoter']=[sum(summary[chrom]['bipromoter'][:i]) for i in range(1,len(summary[chrom]['bipromoter'])+1)]
        summary[chrom]['downstream']=[sum(summary[chrom]['downstream'][:i]) for i in range(1,len(summary[chrom]['downstream'])+1)]

    
    def integrate_summaries(self,summary,chrom):
        """Add the summary of a single chromosome to self.summary['whole']"""
        
        array_adder=lambda x,y: [x[i]+y[i] for i in range(0,len(x))]
        # when self.summary['whole']is empty, just copy the summary of the first chromosome.
        # Then, add array by array

        summary['whole']['promoter']=array_adder(summary['whole']['promoter'],summary[chrom]['promoter'])
        summary['whole']['bipromoter']=array_adder(summary['whole']['bipromoter'],summary[chrom]['bipromoter'])
        summary['whole']['downstream']=array_adder(summary['whole']['downstream'],summary[chrom]['downstream'])
        summary['whole']['gene']=array_adder(summary['whole']['gene'],summary[chrom]['gene'])
        summary['whole']['rel_loc'][0]=array_adder(summary['whole']['rel_loc'][0],summary[chrom]['rel_loc'][0])
        summary['whole']['rel_loc'][1]=array_adder(summary['whole']['rel_loc'][1],summary[chrom]['rel_loc'][1])
        summary['whole']['rel_loc_cds'][0]=array_adder(summary['whole']['rel_loc_cds'][0],summary[chrom]['rel_loc_cds'][0])
        summary['whole']['rel_loc_cds'][1]=array_adder(summary['whole']['rel_loc_cds'][1],summary[chrom]['rel_loc_cds'][1])
        summary['whole']['roi']+=summary[chrom]['roi']
        summary['whole']['Ns']+=summary[chrom]['Ns']

            
    def get_p(self,summary,p,chrom):
        """Get the p of a single chromosome"""
        
        total=summary[chrom]['Ns']
        p.init_table(chrom)
        # check if the denominator is zero
        try:
            p[chrom]['promoter']=map(lambda x: 1.0*x/total,summary[chrom]['promoter'])
        except ZeroDivisionError:
            total=1
            p[chrom]['promoter']=map(lambda x: 1.0*x/total,summary[chrom]['promoter'])
        p[chrom]['bipromoter']=map(lambda x: 1.0*x/total,summary[chrom]['bipromoter'])
        p[chrom]['downstream']=map(lambda x: 1.0*x/total,summary[chrom]['downstream'])
        p[chrom]['gene']=map(lambda x: 1.0*x/total,summary[chrom]['gene'])
            
        #relative locations
        total_rel_loc=sum(summary[chrom]['rel_loc'][0])
        total_rel_loc_cds=sum(summary[chrom]['rel_loc_cds'][0])
        # check if the denominator is zero
        try:
            p_rel_loc0=map(lambda x: 1.0*x/total_rel_loc,summary[chrom]['rel_loc'][0])
        except ZeroDivisionError:
            total_rel_loc=1
            p_rel_loc0=map(lambda x: 1.0*x/total_rel_loc,summary[chrom]['rel_loc'][0])
        p_rel_loc1=map(lambda x: 1.0*x/total_rel_loc,summary[chrom]['rel_loc'][1])
            
        try:
            p_rel_loc_cds0=map(lambda x: 1.0*x/total_rel_loc_cds,summary[chrom]['rel_loc_cds'][0])
        except ZeroDivisionError:
            total_rel_loc_cds=1
            p_rel_loc_cds0=map(lambda x: 1.0*x/total_rel_loc_cds,summary[chrom]['rel_loc_cds'][0])
        p_rel_loc_cds1=map(lambda x: 1.0*x/total_rel_loc_cds,summary[chrom]['rel_loc_cds'][1])
           
        p[chrom]['rel_loc']=[p_rel_loc0,p_rel_loc1]
        p[chrom]['rel_loc_cds']=[p_rel_loc_cds0,p_rel_loc_cds1]
        
        
        try:
            p[chrom]['roi']=1.0*summary[chrom]['roi']/total
        except ZeroDivisionError:
            p[chrom]['roi']=1.0*summary[chrom]['roi']
        except KeyError:
            pass
        
            
        try:
            p[chrom]['chroms']=1.0*summary[chrom]['Ns']/summary['whole']['Ns']
        except ZeroDivisionError:
            p[chrom]['chroms']=total
    
     
    def obtain_distribution_of_sites(self, AnnotT):     
        """Obtain the distribution of sites (eg ChIP regions) over the elements
        
        Parameters:
        1. AnnotT: AnnotTable object (see in tables.py)
        
        Return:
        1. summary: summary of non-overlapping counts of promoter, downstream, gene, and enhancer (and total)
        2. p: proportions of promoter, downstream, gene, and enhancer
        """
        
        # get chromosomes
        chroms = AnnotT.get_chroms()
        if not chroms: return None
        
        summary = {}
        summary['whole'] = {'promoter':[0, 0, 0], 'downstream':[0, 0, 0], 'gene':[0, 0, 0, 0], 'enhancer':0, 'total':0}
        # iterate through chromosomes
        for chrom in chroms:
            summary[chrom]={'promoter':[0, 0, 0], 'downstream':[0, 0, 0], 'gene':[0, 0, 0, 0], 'enhancer':0, 'total':0}
            length = AnnotT.size(chrom)[0]
            for promoter, downstream, gene in itertools.izip(AnnotT[chrom]['promoter'], AnnotT[chrom]['downstream'], AnnotT[chrom]['gene']):
                # get some sort of histograms
                if abs(promoter) == 1:                       # promoter 1
                    summary[chrom]['promoter'][0] += 1
                elif downstream == 1:                   # downstream 1
                    summary[chrom]['downstream'][0] += 1
                elif abs(promoter) == 2:                     # promoter 2
                    summary[chrom]['promoter'][1] += 1
                elif downstream == 2:                   # downstream 2
                    summary[chrom]['downstream'][1] += 1
                elif abs(promoter) == 3:                     # promoter 3
                    summary[chrom]['promoter'][2] += 1
                elif downstream == 3:                   # downstream 3
                    summary[chrom]['downstream'][2] += 1
                elif gene != 0:                         # gene
                    summary[chrom]['gene'][gene-1] += 1
                else:                                       # enhancer
                    summary[chrom]['enhancer'] += 1
            
            # total     
            summary[chrom]['total'] = sum(summary[chrom]['promoter'] + summary[chrom]['downstream'] + summary[chrom]['gene'] + [summary[chrom]['enhancer']])
            
            # update the whole
            summary['whole']['promoter'] = map(lambda x,y: x+y, summary['whole']['promoter'], summary[chrom]['promoter'])
            summary['whole']['downstream'] = map(lambda x,y: x+y, summary['whole']['downstream'], summary[chrom]['downstream'])
            summary['whole']['gene'] = map(lambda x,y: x+y, summary['whole']['gene'], summary[chrom]['gene'])
            summary['whole']['enhancer'] += summary[chrom]['enhancer']
            summary['whole']['total'] += summary[chrom]['total']
    
        p = {}
        chroms += ['whole']
        for chrom in chroms:
            p[chrom] = {'promoter':[0.0, 0.0, 0.0], 'downstream':[0.0, 0.0, 0.0], 'gene':[0.0, 0.0, 0.0, 0.0], 'enhancer':0.0}
            # if total = 0, just go to the next chromosome
            if summary[chrom]['total'] == 0:
                continue
            
            p[chrom]['promoter'] = map(lambda x, y: 1.0*x/y, summary[chrom]['promoter'], [summary[chrom]['total']]*len(summary[chrom]['promoter']))
            p[chrom]['downstream'] = map(lambda x, y: 1.0*x/y, summary[chrom]['downstream'], [summary[chrom]['total']]*len(summary[chrom]['downstream']))
            p[chrom]['gene'] = map(lambda x, y: 1.0*x/y, summary[chrom]['gene'], [summary[chrom]['total']]*len(summary[chrom]['gene']))
            p[chrom]['enhancer'] = 1.0 * summary[chrom]['enhancer'] / summary[chrom]['total']
            
        return summary, p
 
    def obtain_distribution_of_sites_per_chrom(self, AnnotT):     
        """Obtain the distribution of genome  over the elements. 
        Note that this function works only if the annotation has been quantized!!! See def.annotate method.
        
        arguments:
        1. AnnotT: AnnotTable object (see in tables.py)
        
        Return:
        1. summary: summary of non-overlapping counts of promoter, downstream, gene, and enhancer (and total)
        2. p: proportions of promoter, downstream, gene, and enhancer
        """
        
        # get chromosomes
        chroms = AnnotT.get_chroms()
        if not chroms: return None
        
        summary = {}
        #summary['whole'] = {'promoter':[0, 0, 0], 'downstream':[0, 0, 0], 'gene':[0, 0, 0, 0], 'enhancer':0, 'total':0}
        # iterate through chromosomes
        for chrom in chroms:
            summary[chrom]={'promoter':[0, 0, 0], 'downstream':[0, 0, 0], 'gene':[0, 0, 0, 0], 'enhancer':0, 'total':0}
            length = AnnotT.size(chrom)[0]
            for promoter, downstream, gene in itertools.izip(AnnotT[chrom]['promoter'], AnnotT[chrom]['downstream'], AnnotT[chrom]['gene']):
                # get some sort of histograms
                if abs(promoter) == 1:                       # promoter 1
                    summary[chrom]['promoter'][0] += 1
                elif downstream == 1:                   # downstream 1
                    summary[chrom]['downstream'][0] += 1
                elif abs(promoter) == 2:                     # promoter 2
                    summary[chrom]['promoter'][1] += 1
                elif downstream == 2:                   # downstream 2
                    summary[chrom]['downstream'][1] += 1
                elif abs(promoter) == 3:                     # promoter 3
                    summary[chrom]['promoter'][2] += 1
                elif downstream == 3:                   # downstream 3
                    summary[chrom]['downstream'][2] += 1
                elif gene != 0:                         # gene
                    summary[chrom]['gene'][gene-1] += 1
                else:                                       # enhancer
                    summary[chrom]['enhancer'] += 1
            
            # total     
            summary[chrom]['total'] = sum(summary[chrom]['promoter'] + summary[chrom]['downstream'] + summary[chrom]['gene'] + [summary[chrom]['enhancer']])            
        return summary



### test version of annotator for genome background
###
class AnnotatorGBG(Annotator):
    """Annotator inherits Annotator"""
    
    def __init__(self):
        """Constructor"""
        
        Annotator.__init__(self)
        
    def annotate(self,genome_coordinates=None,gene_table=None,roi=None,prom=10000,biprom=20000, down=10000, gene_div=(3,5),quantize=True):
        """Annotate given coordinates based on the given gene table."""
        
        # get the chromsomes of the gene table and genome coordinates
        try:
            chroms_gc=genome_coordinates.keys()
            chroms_gt=gene_table.get_chroms()
            chroms_gt,chroms_gc=set(chroms_gt),set(chroms_gc)
            chroms=chroms_gt.intersection(chroms_gc)
            chroms=list(chroms)
            chroms.sort()
            num_coordinates={}
            num_genes={}
            for chrom in chroms:
                num_coordinates[chrom]=len(genome_coordinates[chrom])
                num_genes[chrom]=len(gene_table[chrom][gene_table[chrom].keys()[0]])
        except AttributeError:
            raise Exception('Genome coordinates and gene table must be given for genome annotation')
        
        
        #initialize with an empty dictionary
        table=AnnotTable()
        #iterate over the chromosomes
        for chrom in chroms: 
            genes=gene_table[chrom]
            num_genes_this_chr=num_genes[chrom]
            coordinates=genome_coordinates[chrom]
            num_coordinates_this_chr=num_coordinates[chrom]
            table.init_table(chrom)
            
            #initialize the promoter distances. This will be used in obtaining bidirectional promoter too.
            prom_dists=[[0,0] for i in xrange(num_coordinates_this_chr)]
        
            # get the nearest genes to set the searching range for promoter and downstream
            nearest_genes=self.find_nearest_genes(genes)
            
            # point begin
            pointerBeg=0
            maxprom = prom
            maxdown = down
            for i in xrange(0,num_genes_this_chr):

                # get the strand of the gene
                try:
                    strand=genes['strand'][i]
                except KeyError:
                    raise Exception("'strand' must be included in the gene annotation table for running CEAS")
                
                # get the beginning and end point of search
                # the beginning and end points are the end of the previous gene and the beginning of the next gene.
                beg,end=0,0
                try:
                    if strand=='+':
                        beg=max(genes['txStart'][i]-maxprom, nearest_genes['before'][i])
                        end=min(genes['txEnd'][i]+maxdown, nearest_genes['after'][i])
                    else:
                        beg=max(genes['txStart'][i]-maxdown, nearest_genes['before'][i])
                        end=min(genes['txEnd'][i]+maxprom, nearest_genes['after'][i])
                except KeyError:    # check the gene annotation table has necessary columns
                    raise Exception("'txStart' and 'txEnd' must be included in the gene annotation table for running CEAS")
                
                # set search index j to the begining point of the last gene. This makes sure that we include isoforms
                j=pointerBeg
                if coordinates[j]>end: continue
                
                # two while loops to detect the annotation start coordinate for the current gene.
                while j>0 and coordinates[j]>=beg: 
                    j-=1  
                while j<num_coordinates_this_chr and coordinates[j]<beg:
                    if j>=table.size(chrom)[0]:
                        table.add_row(chrom,[coordinates[j]]+[0]*table.get_column_num())
                    j+=1
                
                # if get to the end of chromosome, then break
                if j==num_coordinates_this_chr: break
                
                # save the current start point for the next gene
                pointerBeg=j
     
                # otherwise, get the annotations of the probes related with the current gene
                while j<num_coordinates_this_chr and (coordinates[j]>=beg and coordinates[j]<=end):
                    # get the annotation and update the entire annotation table
                    single_annot=self.annotate_single(coordinates[j],strand,genes['txStart'][i],genes['txEnd'][i],\
                                                       genes['cdsStart'][i],genes['cdsEnd'][i],genes['exonStarts'][i],genes['exonEnds'][i],prom,down,gene_div)
                    self.update_annot_table(table,single_annot,coordinates,prom_dists,chrom,j,biprom)
                    j+=1
            
            # quantize promoter, bipromoter and downstream
            if quantize:
                table[chrom]['promoter']=ordinate2(table[chrom]['promoter'],prom)
                table[chrom]['bipromoter']=ordinate2(table[chrom]['bipromoter'],biprom)
                table[chrom]['downstream']=ordinate2(table[chrom]['downstream'],down)
        
        if roi:
            roichroms = roi.get_chroms()
            for chrom in chroms:
                table[chrom]['roi']=[0]*len(genome_coordinates[chrom])
                if chrom in roichroms:
                    self.do_annotation_roi(genome_coordinates[chrom], table[chrom], roi[chrom])
                    
        return table
    
    
    def annotate_single(self,coordinate,strand,txStart,txEnd,cdsStart,cdsEnd,exonStarts,exonEnds,prom,down,gene_div):
        """Annotate a single genome coordinate
        
        Parameters:
        1. coordinate: a list of genome locations to annotate
        2. strand: the strand (+/-) of the gene
        3. txStart: transcription start
        4. txEnd: transcription end
        5. cdsStart: translation start
        6. cdsEnd: translation end
        7. exonStarts: exon start locations
        8. exonEnds: exon end locations
        9. prom: promoter lengths (e.g., (1000, 2000, 3000))
        10. down: downstream lengths (e.g., (1000, 2000, 3000))
        11. gene: the number of divisions of a gene (eg, (3,5))
        
        """
                               
        # container of the annotation for a single location. 
        single_annot=None
        
        maxprom = prom
        maxdown = down
        # get the annotation for the location
        # + strand
        if strand=='+':
            # promoter
            if coordinate<txStart and coordinate>=txStart-maxprom: # txStart-promo <= prob < txStart
                single_annot=['promoter',txStart-coordinate]
            # downstream     
            elif coordinate>=txEnd and coordinate<txEnd+maxdown:# txEnd <= prob < txEnd+down
                single_annot=['downstream',coordinate-txEnd]
            # gene
            elif coordinate>=txStart and coordinate<txEnd:    
                isExon=self.exonify(coordinate,exonStarts,exonEnds) 
                # exon
                if isExon !=None:
                    # coding exon
                    if coordinate>=cdsStart and coordinate<cdsEnd: single_annot=['gene',3,self.get_rel_loc_cds(coordinate,strand,exonStarts,exonEnds,isExon,gene_div)]
                    # 5'UTR
                    elif coordinate<cdsStart: single_annot=['gene',1]
                    # 3'UTR
                    else: single_annot=['gene',2]
                # intron
                else: 
                    single_annot=['gene',4]
                # relative location within the gene
                single_annot.insert(2,self.get_rel_loc(coordinate,strand,txStart,txEnd,gene_div))
        # - strand
        else: 
            # promoter
            if coordinate>=txEnd and coordinate<txEnd+maxprom:  
                single_annot=['promoter',txEnd-coordinate]
            # downstream
            elif coordinate<txStart and coordinate>=txStart-maxdown:
                single_annot=['downstream',txStart-coordinate]
            # gene
            elif coordinate>=txStart and coordinate<txEnd:
                isExon=self.exonify(coordinate,exonStarts,exonEnds) #then in an exon
                # exon
                if isExon != None:
                    # coding exon
                    if coordinate>=cdsStart and coordinate<cdsEnd: single_annot=['gene',3,self.get_rel_loc_cds(coordinate,strand,exonStarts,exonEnds,isExon,gene_div)]
                    # 5'UTR
                    elif coordinate>=cdsEnd: single_annot=['gene',1]
                    # 3'UTR
                    else: single_annot=['gene',2]
                # intron
                else: 
                    single_annot=['gene',4]
                # relative location within the gene
                single_annot.insert(2,self.get_rel_loc(coordinate,strand,txStart,txEnd,gene_div))
        
        return single_annot
    
    
    def update_annot_table(self,table,single_annot,coordinates,prom_dists,chrom,index_coord,biprom):
        """Update the annotation table by comparing the current annotation (single_annot) and the existing annotation in the table"""
        
        categories=['promoter','bipromoter','downstream','gene']
        # If the current annotation is NoneType, which means that the curren genome coordinate does not belong to any of
        # promoters,downstreams,or genes, just add an empty annotation to the end of the annotation table
   
        maxbiprom = biprom
        try:
            category=single_annot[0]
        except TypeError:        # when single_annot=None
            if index_coord>=table.size(chrom)[0]:
                table.add_row(chrom,[coordinates[index_coord]]+[0]*table.get_column_num())
            return
                   
        # When an annotation does not exist at the current genome coordinate, append the new annotation
        # and update
        try:
            original=table[chrom][category][index_coord]
        except IndexError:
                table.add_row(chrom,[coordinates[index_coord]]+[0]*table.get_column_num())
                original=table[chrom][category][index_coord]
    
        # +/- promoter distances
        prom_dist=prom_dists[index_coord]
        
        # when an annotation already exists at the current genome coordinate
#        scaler=lambda dist,maxrange,scale:int(1.0
        new=single_annot[1]
        
        # promoter and bipromoter
        if category=='promoter':
    
            # promoter
            # update + and - distances
            if new < 0:
                if prom_dist[0]==0: 
                    prom_dist[0]=new
                else: 
                    prom_dist[0]=max(new,prom_dist[0])
            elif new > 0:
                if prom_dist[1]==0: 
                    prom_dist[1]=new
                else: 
                    prom_dist[1]=min(new,prom_dist[1])
                        
            # update the promoter annnotation of the current probe
            if original==0:
                if prom_dist[0]==0 or prom_dist[1]==0:
                    fin_prom=abs(sum(prom_dist))
                else: 
                    fin_prom=min(abs(prom_dist[0]),prom_dist[1])
            else: 
                if prom_dist[0]==0 and prom_dist[1]==0:
                    fin_prom=original
                elif prom_dist[0]==0 and prom_dist[1]!=0:
                    fin_prom=min(prom_dist[1],original)
                elif prom_dist[0]!=0 and prom_dist[1]==0:
                    fin_prom=min(abs(prom_dist[0]),original)
                else:
                    fin_prom=min([abs(prom_dist[0]),prom_dist[1],original])
            #update the table with the shortest distance from TSS
            table[chrom][category][index_coord]=fin_prom
            
            # bidirectional promoter
            # if + and - distances have non-zero values, that is bidirectional promoter
            if prom_dist[0]*prom_dist[1] != 0:    
                new_biprom=abs(prom_dist[0])+prom_dist[1]

                #if new_biprom<=cutoff:
                if new_biprom<=maxbiprom:
                    # original bipromoter value
                    ori_biprom=table[chrom]['bipromoter'][index_coord]
                    # update the bipromoter
                    if ori_biprom==0: fin_biprom=new_biprom
                    else: fin_biprom=min(ori_biprom,new_biprom)
                    table[chrom]['bipromoter'][index_coord]=fin_biprom
        elif category=='downstream':
            if original==0: fin_down=new
            else: fin_down=min(new,original)
            table[chrom][category][index_coord]=fin_down
        elif category=='gene':
            if original==0: fin_gene=new
            else: fin_gene=min(new,original)
            table[chrom][category][index_coord]=fin_gene
            table[chrom]['rel_loc'][index_coord]=copy.deepcopy(single_annot[2])
            try:
                table[chrom]['rel_loc_cds'][index_coord]=copy.deepcopy(single_annot[3])
            except IndexError:
                pass
            
                  
    def summarize(self,table,maxprom,maxbiprom,maxdown,by):
        """Provides a summary of the annotation.
        
        Parameters:
        1. table: an AnnotTable object produced by Annotater.annotate()
        2. maxprom: the largest promoter size
        3. maxbiprom: the largest bipromoter size
        4. maxdown: the largest downstream size
        5. by: bin size
        
        Output:
        summary: a Summary object (see tables.py), which contains the summary of the annotation
        """
        
        # get bins from the max sizes of promoters, bidirectional promoters, and downstream
        bins_prom=[0,500]+ range(1000, maxprom, by) + [maxprom]
        bins_biprom=[0,500]+ range(1000, maxbiprom, by) + [maxbiprom]
        bins_down=[0,500]+ range(1000, maxdown, by) + [maxdown]
        
        summary=SummaryGBG(numprom=len(bins_prom)-1,numbiprom=len(bins_biprom)-1,numdown=len(bins_down)-1)
        p=PGBG(numprom=len(bins_prom)-1,numbiprom=len(bins_biprom)-1,numdown=len(bins_down)-1)
          
        try:
            chroms=table.get_chroms()
            if not chroms: raise ValueError
        except (AttributeError,ValueError):
            return summary
               
        # obtain a summary statistics
        summary.init_table('whole')
        for chrom in chroms:
            # get the summary of a single chromosome
            self.get_summary(table,summary,chrom,bins_prom,bins_biprom,bins_down)
            # integrate it into the 'whole'
            self.integrate_summaries(summary,chrom)
            # get p(probablities)
        
        p.init_table('whole')
        for chrom in chroms:
            self.get_p(summary,p,chrom)
        # get the probabilities of the whole genome
        self.get_p(summary,p,'whole')
        
        return summary,p
    
    def get_summary(self,table,summary,chrom,bins_prom,bins_biprom,bins_down):   
        """Get the summary of a single chromosome"""
        
        length=table.size(chrom)[0]
        summary.init_table(chrom)
        summary[chrom]['Ns']=length
                
        # get the binned promoter, bipromoter and downstream distances 
        binned_prom,binned_biprom,binned_down=self.cumbin(table,chrom,bins_prom,bins_biprom,bins_down)
        summary[chrom]['promoter']=binned_prom
        summary[chrom]['bipromoter']=binned_biprom
        summary[chrom]['downstream']=binned_down
                 
        for i in xrange(0,length):
            # get some sort of histograms
             
            if table[chrom]['gene'][i]!=0: summary[chrom]['gene'][table[chrom]['gene'][i]-1]+=1
            if table[chrom]['rel_loc'][i]!=[0,0]: 
                summary[chrom]['rel_loc'][0][table[chrom]['rel_loc'][i][0]-1]+=1
                summary[chrom]['rel_loc'][1][table[chrom]['rel_loc'][i][1]-1]+=1
            if table[chrom]['rel_loc_cds'][i]!=[0,0]:
                summary[chrom]['rel_loc_cds'][0][table[chrom]['rel_loc_cds'][i][0]-1]+=1
                summary[chrom]['rel_loc_cds'][1][table[chrom]['rel_loc_cds'][i][1]-1]+=1
            try:
                if table[chrom]['roi'][i]!=0: summary[chrom]['roi']+=1
            except KeyError:
                pass
        
        # the last gene element is all
        summary[chrom]['gene'][-1]=sum(summary[chrom]['gene'][:-1]) 
        

    def cumbin(self,table,chrom,bins_prom,bins_biprom,bins_down):
        """Calculate the genome bg.
         
        """
    
        allbp = [0] * len(bins_prom)
        allbd = [0] * len(bins_down)
        allbbp = [0] * len(bins_biprom)
    
        prom=table[chrom]['promoter']
        down=table[chrom]['downstream']
        biprom=table[chrom]['bipromoter']
        bp=self.bin(prom, bins_prom)
        bd=self.bin(down, bins_down) 
        bbp=self.bin(biprom, bins_biprom)
        
        ###
        allbp, c = array_adder(allbp, bp)
        allbd, c = array_adder(allbd, bd)
        allbbp, c = array_adder(allbbp, bbp)
        
        # cumsum from 500 to maxprom (or maxbiprom or maxdown)

        cumbp = cumsum(allbp)
        cumbd = cumsum(allbd)
        cumbbp = cumsum(allbbp)
    
        return cumbp, cumbbp, cumbd

    def bin(self, x,bins):
        """Do binning for x.
    
        Parameters:
        1. x: an array of data to do binning for
        2. bins: an array of bins. b[i-1]<= <b[i]
    
        Return:
        1. binned: an array of binned data. Each array element is the count within each bin
        """
    
        #    binlen=len(bins)
        #    binned=[0]*(binlen+1)
        #    sx=sorted(x)
        #    xlen=len(sx)
        #    j=0
        #    
        #    #x < bins[0]
        #    while j<xlen and sx[j]<bins[0]:
        #        binned[0]+=1
        #        j+=1
        #    
        #    # bins[i] <= x < bins[i+1]
        #    for i in xrange(0,binlen-1):
        #        while j<xlen and sx[j]>=bins[i] and sx[j]<bins[i+1]:
        #            binned[i+1]+=1
        #            j+=1
        #    
        #    # x>=bins[-1]
        #    while j < xlen and sx[j]>=bins[-1]:
        #        binned[-1]+=1
        #        j+=1
        #                 
        #    return binned 

        binlen = len(bins)
        binned = [0]* (binlen-1)
        sx = sorted(x)
        xlen = len(sx)
        j=0
    
        # prefilter
        while j < xlen and sx[j]<=bins[0]:
            j+=1
        
        for i in xrange(0, binlen-1):
            while j < xlen and sx[j] > bins[i] and sx[j] <= bins[i+1]:
                binned[i] += 1
                j+=1
    
        return binned


    def obtain_distribution_of_sites_for_genome(self, AnnotT, prom=(1000, 2000, 3000), down=(1000, 2000, 3000)):     
        """Obtain the distribution of genome  over the elements
        
        arguments:
        1. AnnotT: AnnotTable object (see in tables.py)
        2. prom: promoter bins
        3. down: downstream bins
        
        Return:
        1. summary: summary of non-overlapping counts of promoter, downstream, gene, and enhancer (and total)
        2. p: proportions of promoter, downstream, gene, and enhancer
        """
        
        # get chromosomes
        chroms = AnnotT.get_chroms()
        if not chroms: return None
        
        summary = {}
        #summary['whole'] = {'promoter':[0, 0, 0], 'downstream':[0, 0, 0], 'gene':[0, 0, 0, 0], 'enhancer':0, 'total':0}
        # iterate through chromosomes
        for chrom in chroms:
            summary[chrom]={'promoter':[0, 0, 0], 'downstream':[0, 0, 0], 'gene':[0, 0, 0, 0], 'enhancer':0, 'total':0}
            length = AnnotT.size(chrom)[0]
            for promoter, downstream, gene in itertools.izip(AnnotT[chrom]['promoter'], AnnotT[chrom]['downstream'], AnnotT[chrom]['gene']):
                # get some sort of histograms
                if abs(promoter) > 0 and abs(promoter) <= prom[0]:                       # promoter 1
                    summary[chrom]['promoter'][0] += 1
                elif downstream > 0 and downstream <= down[0]:                   # downstream 1
                    summary[chrom]['downstream'][0] += 1
                elif abs(promoter) > prom[0] and abs(promoter) <= prom[1]:                     # promoter 2
                    summary[chrom]['promoter'][1] += 1
                elif downstream > down[0] and downstream <= down[1]:                   # downstream 2
                    summary[chrom]['downstream'][1] += 1
                elif abs(promoter) > prom[1] and abs(promoter) <= prom[2]:                     # promoter 3
                    summary[chrom]['promoter'][2] += 1
                elif downstream > down[1] and downstream <= down[2]:                   # downstream 3
                    summary[chrom]['downstream'][2] += 1
                elif gene != 0:                         # gene
                    summary[chrom]['gene'][gene-1] += 1
                else:                                       # enhancer
                    summary[chrom]['enhancer'] += 1
            
            # total     
            summary[chrom]['total'] = sum(summary[chrom]['promoter'] + summary[chrom]['downstream'] + summary[chrom]['gene'] + [summary[chrom]['enhancer']])            
        return summary

###

class GeneAnnotator:
    """Class GeneAnnotator performs gene-centered annotation given a list of ChIP regions
    
    """
    
    def __init__(self):
        
        self.map = DataFrame()
        self.u = 3000
        self.d = 3000
        
    def annotate(self, GeneT, ChIP, u=3000, d=3000, name2=False):
        """Perform gene-centered annotation
        
        Parameters:
        """
        
        self.u = u
        self.d = d
        
        # initialize the DtaFrame
        self.map.append_column([], colname = 'name')
        self.map.append_column([], colname = 'chr')
        self.map.append_column([], colname = 'txStart')
        self.map.append_column([], colname = 'txEnd')
        self.map.append_column([], colname = 'strand')
        self.map.append_column([], colname = 'dist u TSS')
        self.map.append_column([], colname = 'dist d TSS')
        self.map.append_column([], colname = 'dist u TTS')
        self.map.append_column([], colname = 'dist d TTS')
        self.map.append_column([], colname = '%dbp u TSS' %u)
        self.map.append_column([], colname = '%dbp d TSS' %u)
        self.map.append_column([], colname = '1/3 gene')
        self.map.append_column([], colname = '2/3 gene')
        self.map.append_column([], colname = '3/3 gene')
        self.map.append_column([], colname = '%dbp d TTS' %d)
        self.map.append_column([], colname = 'exons')
        if name2:   # if the user wants to save 'name2'
            self.map.append_column([], colname = 'name2')   
        
        # get the chroms of the gene annotation table
        chroms = GeneT.get_chroms()
        chroms.sort()
        chroms_ChIP = ChIP.get_chroms()
        #maxsearch = 10000   # maximum search range for finding the nearest binding site. 
        
        # iterate through the chromosomes
        for chrom in chroms:
            
            # the number of min distances and the number of annotation fields
            n_mindist = 4           # this variable should be manually updated whenever 'mindists' changes.
            n_annot_fields = 7      # this variable should be manually updated whenever 'annotations' changes.
            
            # if chrom is also in chroms of ChIP, do annotation. Otherwise, just fill out with 'NA'   
            n_genes = len(GeneT[chrom]['txStart'])  
            if chrom in chroms_ChIP:
                n_ChIP = len(ChIP[chrom]['start']) # get the number of genes
                ChIP_start = ChIP[chrom]['start']
                ChIP_end = ChIP[chrom]['end']
                ChIP_center = sorted(map(lambda x, y: (x+y)/2, ChIP_start, ChIP_end))
            else:
                # if this chromosome is not in ChIP, just fill out w/ NA
                mindist = ['NA'] * n_mindist
                annotation = [0] * n_annot_fields
                if name2:
                    for name, chr, txStart, txEnd, strand, n2 in itertools.izip(GeneT[chrom]['name'], [chrom]*n_genes, GeneT[chrom]['txStart'], GeneT[chrom]['txEnd'], GeneT[chrom]['strand'], GeneT[chrom]['name2']):
                        self.map.append_row([name, chr, txStart, txEnd, strand] + mindist + annotation + [n2])
                else:   
                    for name, chr, txStart, txEnd, strand in itertools.izip(GeneT[chrom]['name'], [chrom]*n_genes, GeneT[chrom]['txStart'], GeneT[chrom]['txEnd'], GeneT[chrom]['strand']):
                        self.map.append_row([name, chr, txStart, txEnd, strand] + mindist + annotation)
                
                continue    # then continue to the next chromosome
            
            # 1. get the minimum distances from TSS and TTS            
            # this chromosome's txStart, txEnd, and strand
            txStart = GeneT[chrom]['txStart']
            txEnd = GeneT[chrom]['txEnd']
            strand = GeneT[chrom]['strand']
            
            # get the distances from every txStart to its nearest binding sites in both directions.
            txStart_up, txStart_do = self.min_dists(txStart, ChIP_center)
            
            # get the distances from every txStart to its nearest binding sites in both directions.
            txEnd_up, txEnd_do = self.min_dists(txEnd, ChIP_center)
            
            # re-order the distances of each gene according to the transcription direction
            mindists = self.reorder(txStart_up, txStart_do, txEnd_up, txEnd_do, strand)
                    
            # 2. get the gene-centered annotation
            pointerBeg = 0   
            annotations = []
            
            for txStart, txEnd, strand, exonStarts, exonEnds in itertools.izip(GeneT[chrom]['txStart'], GeneT[chrom]['txEnd'], GeneT[chrom]['strand'], GeneT[chrom]['exonStarts'], GeneT[chrom]['exonEnds']):
                # get the gene region to search for binding site
                if strand == '+':
                    lower = txStart - u 
                    upper = txEnd + d
                else:
                    lower = txStart - d 
                    upper = txEnd + u
                    
                # set search index j to the begining point of the last gene. This makes sure that we include isoforms
                j = pointerBeg
                #if ChIP_start[j] > upper: continue
                
                # adjust the search start point
                while (j > 0 and j < n_ChIP) and ChIP_start[j] >= lower: 
                    j-=1  
                while j < n_ChIP and ChIP_end[j] < lower:
                    j+=1
                
                # if get to the end of ChIP, break
                #if j==n_ChIP: break
                
                # save the current start point for the next gene
                pointerBeg = j
     
                # otherwise, get the annotations of the probes related with the current gene
                while j < n_ChIP and (ChIP_end[j] >= lower and ChIP_start[j] <= upper):
                    j+=1
                        
                # any ChIP region(s) in the search range
                if pointerBeg < j:
                    annotation = self.annotate_single_gene(txStart, txEnd, strand, u, d, exonStarts, exonEnds, ChIP_start[pointerBeg:j], ChIP_end[pointerBeg:j], precision=2)
                else:
                    annotation = [0] * n_annot_fields
                
                annotations.append(annotation)
            
            # 3. save as DataFrame   
            if name2:
                for name, chr, txStart, txEnd, strand, midist, annotation, n2 in itertools.izip(GeneT[chrom]['name'], [chrom]*n_genes, GeneT[chrom]['txStart'], GeneT[chrom]['txEnd'], GeneT[chrom]['strand'], mindists, annotations, GeneT[chrom]['name2']):
                    self.map.append_row([name, chr, txStart, txEnd, strand] + mindist + annotation + [n2])
            else:   
                for name, chr, txStart, txEnd, strand, mindist, annotation in itertools.izip(GeneT[chrom]['name'], [chrom]*n_genes, GeneT[chrom]['txStart'], GeneT[chrom]['txEnd'], GeneT[chrom]['strand'], mindists, annotations):
                    self.map.append_row([name, chr, txStart, txEnd, strand] + mindist + annotation)
            
    
    def min_dists(self, gref, cref):
        """Return the distances to the nearest binding sites upstream and downstream of gene reference point.
        
        ARGUMENTS:
        gref        : a series of gene reference points such as txStart or txEnd
        cref        : a series of ChIP region reference points such as the centers of ChIP regions
            
        """
        
        l_cref = len(cref)
        pointBeg = 0
        up = []
        do = []
        for g in gref:
            #j = bisect.bisect_left(cref, g, pointBeg)
            j = bisect.bisect_left(cref, g)
            
            if j == 0:
                if g == cref[j]:
                    up.append(0)
                    do.append(0)
                else:
                    up.append('NA')
                    do.append(cref[j] - g)
            elif j == l_cref:
                up.append(g - cref[j-1])
                do.append('NA')
            else:
                if g == cref[j]:
                    up.append(0)
                    do.append(0)
                else:
                    up.append(g - cref[j-1])
                    do.append(cref[j] - g)
            
            pointBeg = j
        
        return up, do 
                
                
    def reorder(self, txStart_up, txStart_do, txEnd_up, txEnd_do, strand):
        """Correct the order according to the strand.
        
        ARGUEMNTS:
        txStart_up        : distances to binding sites upstream of txStart
        txStart_do        : distances to binding sites downstream of txStart
        txEnd_up        : distances to binding sites upstream of txEnd
        txEnd_do        : distances to binding sites downstream of txEnd
        strand    : strands of genes
        
        """
        mindists = []
        for txSu, txSd, txEu, txEd, s in itertools.izip(txStart_up, txStart_do, txEnd_up, txEnd_do, strand):
            
            if s == '+':
                mindists.append([txSu, txSd, txEu, txEd])
            else:
                mindists.append([txEd, txEu, txSd, txSu])
        
        return mindists
            
     
    def get_gene_sections(self, txStart, txEnd):
        """Divide equally each gene into 3 sections 
        
        Parameters:
        1. txStart: start position of a gene 
        2. txEnd: end position of the gene
        """
        
        onethird = (txEnd - txStart)/3
        gene = [txStart, txStart + onethird, txEnd - onethird, txEnd]
        
        return gene        
    
    def get_nearTSS_sections(self, txStart, txEnd, strand, u = 3000):
        """Get u bp ustream and downstream of TSS"""
        
        if strand == '+':
            nearTSS = [max(0, txStart - u), txStart, txStart + u]
        else:
            nearTSS = [max(0, txEnd - u), txEnd, txEnd + u]
        
        return nearTSS
    
    def get_nearTTS_sections(self, txStart, txEnd, strand, d = 3000):
        """Get d bp downstream of TTS"""

        if strand == '+':
            nearTTS = [txEnd, txEnd + d]
        else:
            nearTTS = [txStart - d, txStart]
            
        return nearTTS
    
    def extract_txStarts(self,txS,txE,strand,name,sort=True):
        """Extract txStarts given 'txStart', 'txEnd' and 'strand' of a gene annotation table.
        
        Parameters:
        1. txS: 'txStart'
        2. txE: 'txEnd'
        3. starnd: 'strand'
        4. name: 'name'
        4. sort: True=sort by value False=just extract and return
        
        Return:
        a list, refseqs = [(txStart1,strand1,name1),...(txStartn,strandn,namen)] 
        """
        refseqs=[]
        for s,e,st,n in itertools.izip(txS,txE,strand,name):
            if st=='+':
                refseqs.append((s,e,st,n))
            else:
                refseqs.append((e,s,st,n))
        
        if sort:
            refseqs=sorted(refseqs,key=operator.itemgetter(0))
        
        return refseqs      
    
                
    def get_overlap(self, start, end, ChIP_start, ChIP_end):
        """Return the overlap in bp"""
        
        overlap = max(0, min(end, ChIP_end) - max(start, ChIP_start))
        
        return overlap
            
            
    def annotate_single_gene(self, txStart, txEnd, strand, u, d, exonStarts, exonEnds, ChIP_starts, ChIP_ends, precision=None):
        """Annotate the gene with a single ChIP"""
        
        gene = self.get_gene_sections(txStart, txEnd)    # get the sections of a gene to consider
        
        gene_len = len(gene)
        len_in_bp = [gene[i+1] - gene[i] for i in range(gene_len - 1)]
        annot_gene = [0] * (gene_len - 1)
        
        nearTSS = self.get_nearTSS_sections(txStart, txEnd, strand, u)
        nearTSS_len = len(nearTSS)
        len_in_bp_nearTSS = [nearTSS[i+1] - nearTSS[i] for i in range(nearTSS_len - 1)]
        annot_nearTSS = [0] * (nearTSS_len - 1)
        
        nearTTS = self.get_nearTTS_sections(txStart, txEnd, strand, d)
        nearTTS_len = len(nearTTS)
        len_in_bp_nearTTS = [nearTTS[i+1] - nearTTS[i] for i in range(nearTTS_len - 1)]
        annot_nearTTS = [0] * (nearTTS_len - 1)
        
        for ChIP_start, ChIP_end in itertools.izip(ChIP_starts, ChIP_ends):
            temp = [0] * (gene_len - 1)
            for i in range(gene_len - 1):
                temp[i]=self.get_overlap(gene[i], gene[i+1], ChIP_start, ChIP_end)
            annot_gene = map(lambda x, y: x + y, annot_gene, temp)
            
            temp_nearTSS = [0] * (nearTSS_len - 1)
            for i in range(nearTSS_len - 1):
                temp_nearTSS[i] = self.get_overlap(nearTSS[i], nearTSS[i+1], ChIP_start, ChIP_end)
            annot_nearTSS = map(lambda x, y: x + y, annot_nearTSS, temp_nearTSS)
            
            temp_nearTTS = [0] * (nearTTS_len - 1)
            for i in range(nearTTS_len - 1):
                temp_nearTTS[i] = self.get_overlap(nearTTS[i], nearTTS[i+1], ChIP_start, ChIP_end)
            annot_nearTTS = map(lambda x, y: x + y, annot_nearTTS, temp_nearTTS)
            
        # normalize the annotation wrt length 
        annot_gene = [1.0*a/l for a, l in itertools.izip(annot_gene, len_in_bp)]
        annot_nearTSS = [1.0*a/l for a, l in itertools.izip(annot_nearTSS, len_in_bp_nearTSS)]
        annot_nearTTS = [1.0*a/l for a, l in itertools.izip(annot_nearTTS, len_in_bp_nearTTS)]
        
        # if negative strand, reverse the annotation
        if strand =='-':
            annot_gene.reverse()
            annot_nearTSS.reverse()
            annot_nearTTS.reverse()
        
        annot_exons = [0] * len(exonStarts)
        n_exons = len(exonStarts)
        len_exons = sum(map(lambda x, y: y - x, exonStarts, exonEnds))
        # get the coverage on the exons
        for ChIP_start, ChIP_end in itertools.izip(ChIP_starts, ChIP_ends):
            temp_exons = [0] * n_exons
            for i in range(n_exons):
                temp_exons[i] = self.get_overlap(exonStarts[i], exonEnds[i], ChIP_start, ChIP_end)
            annot_exons = map(lambda x, y: x + y, annot_exons, temp_exons)
        
        annot_exons = [1.0 * sum(annot_exons)/ len_exons]
        
        # control the output number precision
        if type(precision) == int:
            if precision > 0:
                annot_nearTSS = map(lambda x: round(x, precision), annot_nearTSS)
                annot_gene = map(lambda x: round(x, precision), annot_gene)
                annot_nearTTS = map(lambda x: round(x, precision), annot_nearTTS)
                annot_exons = map(lambda x: round(x, precision), annot_exons)
            
        return annot_nearTSS + annot_gene + annot_nearTTS + annot_exons
     
    
    def write(self, fn, description=True):
        """Write the gene-centered annotation result in a TXT file with XLS extension
        
        Parameters:
        1. fn: file name. XLS extension will be added automatically.
        2. description: If True, add a header of desciption of columns
        """
        
        # get the span sizes for upstream and downstream
        u = self.u
        d = self.d
        
        # if description is True, put comment (header)
        if description == True:
            
            comment = "\n".join(("# RefSeq: RefSeq ID", \
                                 "# chr: chromosome of a RefSeq gene",\
                                 "# txStart: 5' end of a RefSeq gene", \
                                 "# txEnd: 3' end site of a RefSeq gene", \
                                 "# strand: strand of a RefSeq gene", \
                                 "# dist u TSS: Distance to the nearest ChIP region's center upstream of transcription start site (bp)", \
                                 "# dist d TSS: Distance to the nearest ChIP region's center downstream of transcription start site (bp)",\
                                 "# dist u TTS: Distance to the nearest ChIP region's center upstream of transcription end site (bp)", \
                                 "# dist d TTS: Distance to the nearest ChIP region's center downstream of transcription end (bp)", \
                                 "# %dbp u TSS: Occupancy rate of ChIP region in %dbp upstream of transcription start site (0.0 - 1.0)" %(u, u),\
                                 "# %dbp d TSS: Occupancy rate of ChIP region in %dbp downstream of transcription start site (0.0 - 1.0)" %(u, u), \
                                 "# 1/3 gene: Occupancy rate of ChIP region in 1/3 gene (0.0 - 1.0)", \
                                 "# 2/3 gene: Occupancy rate of ChIP region in 2/3 gene (0.0 - 1.0)", \
                                 "# 3/3 gene: Occupancy rate of ChIP region in 3/3 gene (0.0 - 1.0)", \
                                 "# %dbp d TTS: Occupancy rate of ChIP region in %dbp downstream of transcriptino end (0.0 - 1.0)" %(d, d), \
                                 "# exons: Occupancy rate of ChIP regions in exons (0.0-1.0)", \
                                 "# Note that txStart and txEnd indicate 5' and 3' ends of genes whereas TSS and TTS transcription start and end sites in consideration of strand."))
                                 
            
            comment += "\n"               
        else:
            comment = ""
        
        self.map.write(fn = fn + '.xls', comment = comment)
        


#-------------------------------------
# function
#------------------------------------- 
def make_table_complete(table,chroms):
    """Make the given table complete by adding rows of missing chromosomes.
    
    Some chromosomes might have ChIP regions, which those chromosomes to be empty in the ChIP annotation Summary and P tables.
    In such case, this function fills up the chromosomes with 0 (or 0.0 for P) to prevent errors in advance.
    
    Parameters:
    1. table: a table object (Summary or P, see in tables.py) that will be completed.
    2. chroms: a list of reference chromosomes. Usually from a Summary or P object
    
    """
    
    newchroms = [c for c in chroms if c!= 'whole']
    for chrom in newchroms:
        if not table.has_chrom(chrom): table.init_table(chrom)
        
    
def estimate_pvals(genome_p,ChIP_summary,ChIP_p):
    """Estimate p values using the binomial modeling
    
    Parameters:
    1. genome_p: a P object (see tables.py) of genome background statistics.
    2. ChIP_summary: a Summary object (see tables.py), which contains the summary of ChIP annotation.
    3. ChIP_p: a P object (see tables.py), which contains the probabilities of ChIP annotation.
    Output:
    pvalue: a P object of p-values.
    """
    
    pvalue=P()
    chroms=set(genome_p.get_chroms()).intersection(set(ChIP_summary.get_chroms()))
    chroms = list(chroms)
    
    
    for chrom in chroms:    
        _get_pvals(genome_p,ChIP_summary,ChIP_p,pvalue,chrom)
    
    N=ChIP_summary['whole']['Ns']
    for chrom in chroms:
        q=ChIP_summary[chrom]['Ns']
        p=genome_p[chrom]['chroms']
        pa=ChIP_p[chrom]['chroms']
        pvalue[chrom]['chroms']=_calc_pval(q,N,p,pa)
        
    return pvalue
        
def _get_pvals(genome_p,ChIP_summary,ChIP_p,pvalue,chrom):   
    """Get pvalues for the given chromosome"""
    
    pvalue.init_table(chrom)
    length=len(ChIP_summary[chrom]['promoter'])
    N=ChIP_summary[chrom]['Ns']
    for i in range(0,length):
        q=ChIP_summary[chrom]['promoter'][i]
        p=genome_p[chrom]['promoter'][i]
        pa=ChIP_p[chrom]['promoter'][i]
        pvalue[chrom]['promoter'][i]=_calc_pval(q,N,p,pa)
    
    length=len(ChIP_summary[chrom]['bipromoter'])
    N=ChIP_summary[chrom]['Ns']
    for i in range(0,length):
        q=ChIP_summary[chrom]['bipromoter'][i]
        p=genome_p[chrom]['bipromoter'][i]
        pa=ChIP_p[chrom]['bipromoter'][i]
        pvalue[chrom]['bipromoter'][i]=_calc_pval(q,N,p,pa)
        
    length=len(ChIP_summary[chrom]['downstream'])
    N=ChIP_summary[chrom]['Ns']
    for i in range(0,length):
        q=ChIP_summary[chrom]['downstream'][i]
        p=genome_p[chrom]['downstream'][i]
        pa=ChIP_p[chrom]['downstream'][i]
        pvalue[chrom]['downstream'][i]=_calc_pval(q,N,p,pa)
    
    length=len(ChIP_summary[chrom]['gene'])
    N=ChIP_summary[chrom]['Ns']
    for i in range(0,length):
        q=ChIP_summary[chrom]['gene'][i]
        p=genome_p[chrom]['gene'][i]
        pa=ChIP_p[chrom]['gene'][i]
        pvalue[chrom]['gene'][i]=_calc_pval(q,N,p,pa)
    
    length=len(ChIP_summary[chrom]['rel_loc'][0])
    N=sum(ChIP_summary[chrom]['rel_loc'][0])
    for i in range(0,length):
        q=ChIP_summary[chrom]['rel_loc'][0][i]
        p=genome_p[chrom]['rel_loc'][0][i]
        pa=ChIP_p[chrom]['rel_loc'][0][i]
        pvalue[chrom]['rel_loc'][0][i]=_calc_pval(q,N,p,pa)
    
    length=len(ChIP_summary[chrom]['rel_loc'][1])
    for i in range(0,length):
        q=ChIP_summary[chrom]['rel_loc'][1][i]
        p=genome_p[chrom]['rel_loc'][1][i]
        pa=ChIP_p[chrom]['rel_loc'][1][i]
        pvalue[chrom]['rel_loc'][1][i]=_calc_pval(q,N,p,pa)

    length=len(ChIP_summary[chrom]['rel_loc_cds'][0])
    N=sum(ChIP_summary[chrom]['rel_loc_cds'][0])
    for i in range(0,length):
        q=ChIP_summary[chrom]['rel_loc_cds'][0][i]
        p=genome_p[chrom]['rel_loc_cds'][0][i]
        pa=ChIP_p[chrom]['rel_loc_cds'][0][i]
        pvalue[chrom]['rel_loc_cds'][0][i]=_calc_pval(q,N,p,pa)
    
    length=len(ChIP_summary[chrom]['rel_loc_cds'][1])
    for i in range(0,length):
        q=ChIP_summary[chrom]['rel_loc_cds'][1][i]
        p=genome_p[chrom]['rel_loc_cds'][1][i]
        pa=ChIP_p[chrom]['rel_loc_cds'][1][i]
        pvalue[chrom]['rel_loc_cds'][1][i]=_calc_pval(q,N,p,pa)
    

    N=ChIP_summary[chrom]['Ns']
    q=ChIP_summary[chrom]['roi']
    p=genome_p[chrom]['roi']
    pa=ChIP_p[chrom]['roi']
    pvalue[chrom]['roi']=_calc_pval(q,N,p,pa)
    
    N=ChIP_summary['whole']['Ns']
    q=ChIP_summary[chrom]['Ns']
    p=genome_p[chrom]['chroms']
    pa=ChIP_p[chrom]['chroms']
    pvalue[chrom]['chroms']=_calc_pval(q,N,p,pa)
            
def _calc_pval(q,N,p,pa):
    """Calculate a pvalue given N,q,p,pa"""
    
    if p>=pa:
        pval=Prob.binomial_cdf(q,N,p)
    else: 
        pval=Prob.binomial_cdf(q,N,p,lower=False)
        
    if pval==0.0: pval=4.92e-324
    
    return pval

def estimate_enhancer_p(summary):
    """Estimate the proportion of intergenic enhancer
    
    Intergenic enhancers are defined as the remaining part after removing promoter, downstream, and genic regions
    """

    
