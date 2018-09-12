#!/usr/bin/env python3

# various functions and mixins for downstream genomic and epigenomic anlyses

import os
import glob
import re
import random
from datetime import datetime
import time
import subprocess
from math import floor,log10

from pybedtools import BedTool
import pandas as pd
import gseapy
import rpy2.robjects as ro
import rpy2.rinterface as ri
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from IPython.display import Image,display

#Get Current Git Commit Hash for version
path = [x.replace(' ','\ ') for x in os.popen('echo $PYTHONPATH').read().split(':') if 'dkfunctions' in x.split('/')]

if len(path) > 0:
    version = os.popen('cd {}; git rev-parse HEAD'.format(path[0])).read()[:-1]
    __version__ = 'v0.1, Git SHA1: {}'.format(version)
else:
    __version__ = 'v0.1, {:%Y-%m-%d}'.format(datetime.now())

def read_pd(file,*args,**kwargs):
    if (file.split('.')[-1] == 'txt') or (file.split('.')[-1] == 'tab'):
        return pd.read_table(file, header= 0, index_col=0, *args, **kwargs)
    elif (file.split('.')[-1] == 'xls') or (file.split('.')[-1] == 'xlsx'):
        return pd.read_excel(file, *args, **kwargs)
    else:
        raise IOError("Cannot parse count matrix.  Make sure it is .txt, .xls, or .xlsx")

def rout_write(x):
    '''
    function for setting r_out to print to file instead of jupyter
    rpy2.rinterface.set_writeconsole_regular(rout_write)
    rpy2.rinterface.set_writeconsole_warnerror(rout_write)
    '''
    print(x, file=open('{}/R_out_{:%Y-%m-%d}.txt'.format(os.getcwd(), datetime.now()), 'a'))

def annotate_peaks(dict_of_dfs, folder, genome, db='UCSC', check=False):
    '''
    Annotate a dictionary of dataframes from bed files to the genome using ChIPseeker and Ensembl annotations.

    Inputs
    ------
    dict_of_beds: dictionary of bed files
    folder: output folder
    genome: hg38, hg19, mm10
    db: default UCSC, but can also accept Ensembl
    check: bool. checks whether annotation file already exists

    Returns
    -------
    dictionary of annotated bed files as dataframe

    '''
    pandas2ri.activate()

    ri.set_writeconsole_regular(rout_write)
    ri.set_writeconsole_warnerror(rout_write)

    chipseeker = importr('ChIPseeker')
    genomicFeatures = importr('GenomicFeatures')
    makeGR = ro.r("makeGRangesFromDataFrame")
    as_df = ro.r("as.data.frame")    

    check_df = {key:os.path.isfile('{}{}_annotated.txt'.format(folder,key.replace(' ','_'))) for key in dict_of_dfs.keys()}
    return_bool = False not in set(check_df.values())
    if return_bool & check:
        return {'{}_annotated'.format(key):pd.from_csv('{}{}_annotated.txt'.format(folder,key.replace(' ','_')), index_col=0, header=0, sep="\t") for key in dict_of_dfs.keys()}

    if db.lower() == 'ucsc':
        species = ('Mmusculus' if genome.lower() == 'mm10' else 'Hsapiens')
        TxDb = importr('TxDb.{}.UCSC.{}.knownGene'.format(species, genome.lower()))
        txdb = ro.r('txdb <- TxDb.{}.UCSC.{}.knownGene'.format(species, genome.lower()))
    elif db.lower() == 'ensembl':
        pwd = '/Volumes/ce-bioinformatics/COMMON/txdb/gencode_{}_txdb_2018_05_23.sqlite'
        loadDb = ro.r('loadDb')
        txdb = loadDb(pwd.format(genome.lower()))
    else:
        raise ValueError('UCSC or Ensembl only.')

    os.makedirs(folder, exist_ok=True)
    
    if genome.lower() == 'mm10':
        annoDb = importr('org.Mm.eg.db')
        anno = 'org.Mm.eg.db'
    elif genome.lower() == 'hg38' or genome.lower() == 'hg19':
        annoDb = importr('org.Hs.eg.db')
        anno = 'org.Hs.eg.db'

    return_dict={}
    
    print('Annotating Peaks...')
    for key, df in dict_of_dfs.items():
        if check & check_df[key]:
            return_dict['{}_annotated'.format(key)] = pd.from_csv('{}{}_annotated.txt'.format(folder,key.replace(' ','_')), index_col=0, header=0, sep="\t")
        else:
            col_len=len(df.columns)
            df.columns = ["chr","start","end"] + list(range(col_len - 3))
            GR = makeGR(df)
            GR_anno = chipseeker.annotatePeak(GR, overlap='all', TxDb=txdb, annoDb=anno)
            return_dict['{}_annotated'.format(key)] = ro.pandas2ri.ri2py(chipseeker.as_data_frame_csAnno(GR_anno))
            return_dict['{}_annotated'.format(key)].to_csv('{}{}_annotated.txt'.format(folder,key.replace(' ','_')),index=True, header=True, sep="\t")
    
    return return_dict

def plot_peak_genomic_annotation(dict_of_df,folder,genome):
    '''
    Plots genomic annotation graphs.

    Inputs
    ------
    dict_of_df: dictionary of dataframes of bed type data
    folder: output folder

    Returns
    -------
    NoneType

    '''
    pandas2ri.activate()
    ri.set_writeconsole_regular(rout_write)
    ri.set_writeconsole_warnerror(rout_write)
    
    species = ('Mmusculus' if genome.lower() == 'mm10' else 'Hsapiens')

    TxDb = importr('TxDb.{}.UCSC.{}.knownGene'.format(species, genome.lower()))
    txdb = ro.r('txdb <- TxDb.{}.UCSC.{}.knownGene'.format(species, genome.lower()))

    if genome.lower() == 'mm10':
        annoDb = importr('org.Mm.eg.db')
        anno = 'org.Mm.eg.db'
    elif genome.lower() == 'hg38' or genome.lower() == 'hg19':
        annoDb = importr('org.Hs.eg.db')
        anno = 'org.Hs.eg.db'

    chipseeker = ro.packages.importr('ChIPseeker')
    grdevices = ro.packages.importr('grDevices')
    annotatePeak = ro.r('annotatePeak')
    makeGR = ro.r("makeGRangesFromDataFrame")
    as_df = ro.r("as.data.frame")
    c = ro.r('c')
    lapply= ro.r('lapply')
    upsetplot = ro.r('upsetplot')
    plot = ro.r('plot')
    
    out = folder + 'all_peak_annotation'
    os.makedirs(out, exist_ok=True)
    
    GR={}
    GR_anno={}
    for key, df in dict_of_df.items():
        col_len=len(df.columns)
        df.columns = ["chr","start","end"] + list(range(col_len - 3))
        GR[key] = makeGR(df)
        GR_anno[key] = chipseeker.annotatePeak(GR[key], overlap='all', TxDb=txdb, annoDb="org.Hs.eg.db")

        grdevices.png(file='{}/{}_annoBar.png'.format(out,key.replace(' ','_')), width=512, height=256)
        plot(chipseeker.plotAnnoBar(GR_anno[key]))
        grdevices.dev_off()

        grdevices.png(file='{}/{}_TSS_Bar.png'.format(out,key.replace(' ','_')), width=512, height=256)
        plot(chipseeker.plotDistToTSS(GR_anno[key]))
        grdevices.dev_off()

        grdevices.png(file='{}/{}_annoData.png'.format(out,key).replace(' ','_'), width=1000, height=500)
        upsetplot(GR_anno[key], vennpie=True)
        grdevices.dev_off()
    
def plot_venn2(Series, string_name_of_overlap, folder):
    '''
    Series with with overlaps 10,01,11
    Plots a 2 way venn.
    Saves to file.
    '''
    
    folder = '{}venn2/'.format(folder) if folder.endswith('/') else '{}/venn2/'.format(folder)
    os.makedirs(folder, exist_ok=True)
    
    plt.figure(figsize=(7,7))
    
    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
           }
    
    plt.rc('font', **font)
  
    #make venn
    venn_plot = venn2(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]), set_labels = Series.index.tolist())
    patch=['10','01','11']
    colors=['green','blue','teal']
    for patch,color in zip(patch,colors):
        venn_plot.get_patch_by_id(patch).set_color('none')
        venn_plot.get_patch_by_id(patch).set_alpha(.4)
        venn_plot.get_patch_by_id(patch).set_edgecolor('none')
    
    c= venn2_circles(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]))
    colors_test=['green','blue']
    for circle,color in zip(c,colors_test):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(3)
     
    plt.title(string_name_of_overlap + " overlaps")
    plt.tight_layout()
    plt.savefig('{}{}-overlap.svg'.format(folder, string_name_of_overlap.replace(' ','_')))
    plt.savefig('{}{}-overlap.png'.format(folder, string_name_of_overlap.replace(' ','_')), dpi=300)
    
def plot_venn2_set(dict_of_sets, string_name_of_overlap, folder):
    '''
    Plots a 2 way venn from a dictionary of sets
    Saves to file.

    Inputs
    ------
    dict_of_sets: dictionary of sets to overlap
    string_name_of_overlap: string with name of overlap
    folder: output folder

    Returns
    -------
    None

    '''
    folder = '{}venn2/'.format(folder) if folder.endswith('/') else '{}/venn2/'.format(folder)
    os.makedirs(folder, exist_ok=True)
    
    plt.figure(figsize=(7,7))
    
    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
           }
    
    plt.rc('font', **font)
  
    set_list = []
    set_names = []
    for name,setlist in dict_of_sets.items():
        set_list.append(setlist)
        set_names.append(name)
        
    #make venn
    venn_plot = venn2(subsets=set_list, set_labels = set_names)
    patch=['10','01','11']
    colors=['green','blue','teal']
    for patch,color in zip(patch,colors):
        venn_plot.get_patch_by_id(patch).set_color('none')
        venn_plot.get_patch_by_id(patch).set_alpha(.4)
        venn_plot.get_patch_by_id(patch).set_edgecolor('none')
    
    c= venn2_circles(subsets=set_list)
    colors_test=['green','blue']
    for circle,color in zip(c,colors_test):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(3)
     
    plt.title(string_name_of_overlap + " overlaps")
    plt.tight_layout()
    plt.savefig('{}{}-overlap.svg'.format(folder, string_name_of_overlap.replace(' ','_')))
    plt.savefig('{}{}-overlap.png'.format(folder, string_name_of_overlap.replace(' ','_')), dpi=300)

def plot_venn3_set(dict_of_sets, string_name_of_overlap, folder):
    '''  
    Makes 3 way venn from 3 sets.
    Saves to file.

    Inputs
    ------
    dict_of_sets: dictionary of sets to overlap
    string_name_of_overlap: string with name of overlap
    folder: output folder

    Returns
    -------
    None

    '''
    folder = '{}venn3/'.format(folder) if folder.endswith('/') else '{}/venn3/'.format(folder)
    os.makedirs(folder, exist_ok=True)

    plt.figure(figsize=(7,7))
    
    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
           }
    
    plt.rc('font', **font)
    
    set_list = []
    set_names = []
    for name,setlist in dict_of_sets.items():
        set_list.append(setlist)
        set_names.append(name)
    
    #make venn
    venn_plot = venn3(subsets=set_list, set_labels = set_names)
    patch=['100','110','101','010','011','001','111']
    for p in patch:
        if venn_plot.get_patch_by_id(p):
            venn_plot.get_patch_by_id(p).set_color('none')
            venn_plot.get_patch_by_id(p).set_alpha(.4)
            venn_plot.get_patch_by_id(p).set_edgecolor('none')
    
    #make
    c= venn3_circles(subsets=set_list)
    colors_list=['green','blue','grey']
    for circle,color in zip(c,colors_list):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(4)
     
    plt.title("{} Overlaps".format(string_name_of_overlap))
    plt.tight_layout()
    plt.savefig('{}{}-overlap.svg'.format(folder, string_name_of_overlap.replace(' ','_')))
    plt.savefig('{}{}-overlap.png'.format(folder, string_name_of_overlap.replace(' ','_')), dpi=300)
    
def plot_venn3_counts(element_list, set_labels, string_name_of_overlap, folder):
    '''    
    Plot three way venn based on counts of specific overlaping numbers.
    Saves to file.

    Inputs
    ------
    element_list: tuple with counts of the the overlaps from (Abc,aBc,ABc,abC,AbC,ABC)
    set_labels: list or tuple with names of the overlaps ('A','B','C')
    string_name_of_overlap: string with name of overlap
    folder: output folder

    Returns
    -------
    None

    '''
    folder = '{}venn3/'.format(folder) if folder.endswith('/') else '{}/venn3/'.format(folder)
    os.makedirs(folder, exist_ok=True)
    
    
    plt.figure(figsize=(7,7))
    
    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
           }
    
    plt.rc('font', **font)
    
    #make venn
    venn_plot = venn3(subsets=element_list, set_labels = set_labels)
    patch=['100','110','101','010','011','001','111']
    for p in patch:
        if venn_plot.get_patch_by_id(p):
            venn_plot.get_patch_by_id(p).set_color('none')
            venn_plot.get_patch_by_id(p).set_alpha(.4)
            venn_plot.get_patch_by_id(p).set_edgecolor('none')
    
    #make
    c= venn3_circles(subsets=element_list)
    colors_list=['green','blue','grey']
    for circle,color in zip(c,colors_list):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(4)
     
    plt.title("{} Overlaps".format(string_name_of_overlap))
    plt.tight_layout()
    plt.savefig('{}{}-overlap.svg'.format(folder, string_name_of_overlap.replace(' ','_')))
    plt.savefig('{}{}-overlap.png'.format(folder, string_name_of_overlap.replace(' ','_')), dpi=300)
    
def overlap_two(bed_dict, genome=None):
    '''
    Takes a dictionary of two bed-like format files.
    Merges all overlapping peaks for each bed into a master file.
    Intersects beds to merged master file.
    Performs annotations with ChIPseeker if genome is specified.
    Plots venn diagrams of peak overlaps
    If genome is specified, also plots venn diagrams of annotated gene sets.

    Inputs
    ------
    bed_dict:  dictionary of BedTool files
    genome: 'hg38','hg19','mm10'
    
    Returns
    -------
    Returns a dictionary of dataframes from unique and overlap peaks.
    If genome is specified, includes a dictionary of annotated peaks.
    '''

    names = list(bed_dict.keys())

    Folder = '{}/'.format(os.getcwd())
    subfolder = '{}_{}_overlap/'.format(names[0].replace(' ','_'), names[1].replace(' ','_'))
    
    out = '{}{}'.format(Folder, subfolder)
    os.makedirs(out, exist_ok=True)
    print('Output files are found in {}'.format(out))
    
    masterfile = bed_dict[names[0]].cat(bed_dict[names[1]]).sort().merge()
    sorted_dict = {key:bed.sort().merge() for key, bed in bed_dict.items()}
    overlap_dict = {'overlap':masterfile.intersect(sorted_dict[names[0]]).intersect(sorted_dict[names[1]])}
    for key,bed in sorted_dict.items():
        other = {other_key:other_bed for other_key,other_bed in sorted_dict.items() if other_key != key}
        overlap_dict[key] = masterfile.intersect(sorted_dict[key]).intersect(list(other.values())[0], v=True)

    for key,bed in overlap_dict.items():
        bed.to_dataframe().to_csv('{}{}{}-unique-peaks-from-mergedPeaks.bed'.format(Folder, subfolder, key.replace(' ','_')),
                                  header=None, index=None, sep="\t")
    
    overlap_numbers = pd.Series({names[0]: len(overlap_dict[names[0]]),
                                 names[1]: len(overlap_dict[names[1]]),
                                 'overlap': len(overlap_dict['overlap'])
                                 },
                                index=[names[0], names[1],'overlap']         
                                )

    #Venn
    plot_venn2(overlap_numbers, 
               '{} and\n{} peak'.format(names[0], names[1]),
               '{}{}'.format(Folder,subfolder)
              )
    if bool(genome):
        print('Annotating overlaping peaks...')
        #Annotate with ChIPseeker
        unikey='{}_unique'
        unianno='{}_unique_annotated'
        return_dict = annotate_peaks({unikey.format(key):bed.to_dataframe() for key,bed in overlap_dict.items()}, '{}{}'.format(Folder,subfolder), genome=genome)

        Set1_unique = set(return_dict[unianno.format(names[0])].SYMBOL.unique().tolist())
        Set2_unique = set(return_dict[unianno.format(names[1])].SYMBOL.unique().tolist()) 
        Overlap_Set = set(return_dict[unianno.format('overlap')].SYMBOL.unique().tolist())

        venn2_dict = {names[0]: (Set1_unique | Overlap_Set),
                      names[1]: (Set2_unique | Overlap_Set)
                     }

        plot_venn2_set(venn2_dict,
                      '{} and {}\nannotated gene'.format(names[0], names[1]),
                      '{}{}'.format(Folder,subfolder)
                      )
        
        gene_overlaps = {}
        gene_overlaps['{}_unique_genes'.format(names[0])] = Set1_unique - (Set2_unique | Overlap_Set)
        gene_overlaps['{}_unique_genes'.format(names[1])] = Set2_unique - (Set1_unique | Overlap_Set)
        gene_overlaps['Overlap_Gene_Set'] = (Set1_unique & Set2_unique) | Overlap_Set
        
        #for key, gene_set in gene_overlaps.items():
        #    out_file = '{}{}{}.txt'.format(Folder,subfolder,key)
        #    for gene in gene_set:
        #        print(gene, file=open(out_file, 'a'))
        
        for key,item in gene_overlaps.items():
            return_dict[key]=item
            
        for key,df in overlap_dict.items():
            return_dict[key]=df
        
    else: 
        return_dict = overlap_dict
    
    return return_dict

def enrichr(gene_list, description, out_dir):
    '''
    Performs GO Molecular Function, GO Biological Process and KEGG enrichment on a gene list.
    Uses enrichr.

    Inputs
    ------
    gene_list: list of genes to perform enrichment on
    description: string description for title
    out_dir: output director

    Returns
    -------

    None

    '''    
    os.makedirs(out_dir, exist_ok=True)
    
    gseapy.enrichr(gene_list=gene_list,
                   description='{}_KEGG'.format(description),
                   gene_sets='KEGG_2016', 
                   outdir=out_dir
                   )
    gseapy.enrichr(gene_list=gene_list,
                   description='{}_GO_biological_process'.format(description),
                   gene_sets='GO_Biological_Process_2017b', 
                   outdir=out_dir
                  )
    gseapy.enrichr(gene_list=gene_list, 
                   description='{}_GO_molecular_function'.format(description),
                   gene_sets='GO_Molecular_Function_2017b', 
                   outdir=out_dir
                  )

def overlap_three(bed_dict, genome=None):
    '''
    Takes a dictionary of three bed-like format files.
    Merges all overlapping peaks for each bed into a master file.
    Intersects beds to merged master file.
    Performs annotations with ChIPseeker if genome is specified.
    Plots venn diagrams of peak overlaps
    If genome is specified, also plots venn diagrams of annotated gene sets.

    Inputs
    ------
    bed_dict:  dictionary of BedTool files
    genome: 'hg38','hg19','mm10'
    
    Returns
    -------
    Returns a dictionary of dataframes from unique and overlap peaks.
    If genome is specified, includes a dictionary of annotated peaks.
    '''
    from collections import OrderedDict

    names = list(bed_dict.keys())

    Folder = '{}/'.format(os.getcwd())
    subfolder = '{}-{}-{}-overlap/'.format(names[0].replace(' ','_'), names[1].replace(' ','_'), names[2].replace(' ','_'))
    
    out = '{}{}'.format(Folder, subfolder)
    os.makedirs(out, exist_ok=True)
    print('Output files are found in {}'.format(out))
    print('A: {}, B: {}, C: {}'.format(names[0], names[1], names[2]))
    
    master = bed_dict[names[0]].cat(bed_dict[names[1]]).cat(bed_dict[names[2]]).sort().merge()
    
    A = bed_dict[names[0]].sort().merge()
    B = bed_dict[names[1]].sort().merge()
    C = bed_dict[names[2]].sort().merge()

    sorted_dict=OrderedDict({'master':master,'A':A,'B':B,'C':C})
    sorted_dict['Abc']=master.intersect(A).intersect(B, v=True).intersect(C,v=True)
    sorted_dict['aBc']=master.intersect(B).intersect(A, v=True).intersect(C,v=True)
    sorted_dict['ABc']=master.intersect(A).intersect(B).intersect(C,v=True)
    sorted_dict['abC']=master.intersect(C).intersect(A, v=True).intersect(B,v=True)
    sorted_dict['AbC']=master.intersect(A).intersect(C).intersect(B, v=True)
    sorted_dict['aBC']=master.intersect(B).intersect(C).intersect(A,v=True)
    sorted_dict['ABC']=master.intersect(A).intersect(B).intersect(C)

    labTup =tuple(key for key in sorted_dict.keys())
    lenTup =tuple(len(bed) for bed in sorted_dict.values())
    
    print('{}\n{}'.format(labTup, lenTup))

    plot_venn3_counts(lenTup[4:], names, '', out)

    for key,bed in sorted_dict.items():
        bed.to_dataframe().to_csv('{}{}-peaks-from-mergedPeaks.bed'.format(out, key.replace(' ','_')),
                                  header=None, index=None, sep="\t")

    if bool(genome):
        print('Annotating ovelapped peaks...')
        unikey='{}_unique'
        unianno='{}_unique_annotated'
        return_dict = annotate_peaks({unikey.format(key):bed.to_dataframe() for key,bed in sorted_dict.items()}, out, genome=genome)

        Set1=set(return_dict[unianno.format('A')].SYMBOL.unique().tolist())
        Set2=set(return_dict[unianno.format('B')].SYMBOL.unique().tolist())
        Set3=set(return_dict[unianno.format('C')].SYMBOL.unique().tolist())

        plot_venn3_set({names[0]:Set1, names[1]:Set2, names[2]:Set3}, '{}_{}_{}'.format(names[0],names[1],names[2]), out)

    return sorted_dict if genome == None else {**sorted_dict, **return_dict}

def splice_bar(data, title, x, y):
    
    '''
    Plots bar graph of misplicing counts as file.

    Inputs
    ------
    data: dataframe
    title: string plot title
    x: string of columm title for number of events in data
    y: string of column title for splicing type in data

    Returns
    -------
    None
    '''
    sns.set(context='paper', font='Arial', style='white', font_scale=2)

    plot = sns.barplot(x = x, y=y, data = data)
    plot.set_title(title)
    plot.set_ylabel('')

    sns.despine()
    sns.utils.plt.savefig('{}.png'.format(title.replace(' ','_')), dpi=300)

def make_df(dict_of_sets, name):
    '''
    Make a dataframe from a dictionary of sets.

    Inputs
    ------
    dict_of_sets: dictionary of sets
    name: string name of file

    Returns
    -------

    dataframe

    '''

    out_dir='{pwd}/{name}/'.format(pwd=os.getcwd(), name=name.replace(' ','_'))
    os.makedirs(out_dir, exist_ok=True)
    
    count = 0
    for key,genes in dict_of_sets.items():
        count = max(count,len(genes))
    
    df = pd.DataFrame(index = range(1,count+1))

    for key,genes in dict_of_sets.items():
        df[key] = pd.Series(list(genes) + ['NA']*(count-len(genes)))
    
    df.to_excel('{}/{}.xls'.format(out_dir,name.replace(' ','_')), index=False)

    return df

def enrichr_topterm(gene_list, description, out_dir, top_term, figsize):
    '''
    Performs GO Molecular Function, GO Biological Process and KEGG enrichment on a gene list.
    Uses enrichr.
    plots the top specified number of terms

    Inputs
    ------
    gene_list: list of genes to perform enrichment on
    description: string description for title
    out_dir: output directory
    top_term: integer, number of top hits to plot
    figsize: tuple, figure size.  (6,12)

    Returns
    -------

    None

    '''
    
    os.makedirs(out_dir, exist_ok=True)
    
    gseapy.enrichr(figsize= figsize,
                   top_term=top_term,
                   gene_list=gene_list,
                   description='{}_KEGG'.format(description),
                   gene_sets='KEGG_2016', 
                   outdir=out_dir
                   )

    gseapy.enrichr(figsize= figsize,
                   top_term=top_term,
                   gene_list=gene_list,
                   description='{}_GO_biological_process'.format(description),
                   gene_sets='GO_Biological_Process_2017b', 
                   outdir=out_dir
                  )
    gseapy.enrichr(figsize= figsize,
                   top_term=top_term,
                   gene_list=gene_list, 
                   description='{}_GO_molecular_function'.format(description),
                   gene_sets='GO_Molecular_Function_2017b', 
                   outdir=out_dir
                  )

def plot_col(df, title, ylabel, out='', xy=(None,None), xticks=None, plot_type=['violin'], pvalue=False, compare_tags=None):
    '''
    Two column boxplot from dataframe.  Titles x axis based on column names.
    
    Inputs
    ------
    df: dataframe (uses first two columns)
    title: string of title
    ylabel: string of y label
    xy: If specified, will x is the label column and y is the data column. (default: (None,None): Data separated into two columns).
    xticks: list of xtick names (default is column name)
    pvalue: bool to perform ttest (default is False).  Will only work if xy=(None,None) or ther are only two labels in x. 
    plot_type: list of one or more: violin, box, swarm (default=violin)
    compare_tags:  if xy and pvalue is specified and there are more than two tags in x, specify the tags to compare. eg. ['a','b']
    out: out parent directory.  if none returns into colplot/

    Returns
    ------
    None
    
    
    '''
    if len(out) != 0:
        out = out if out.endswith('/') else '{}/'.format(out)
    out = '{}colplot/'.format(out)
    os.makedirs(out, exist_ok=True)

    plt.clf()
    sns.set(context='paper', font='Arial', font_scale=2, style='white', rc={'figure.dpi': 300, 'figure.figsize':(5,6)})
    
    if type(plot_type) != list:
        plot_type = plot_type.split()
    lower_plot_type = [x.lower() for x in plot_type]

    if len(lower_plot_type) == 0:
    	raise IOError('Input a plot type.')
    elif True not in {x in lower_plot_type for x in ['violin', 'box', 'swarm']}:
        raise IOError('Did not recognize plot type.')

    if 'swarm' in lower_plot_type:
        if xy == (None,None):
            fig = sns.swarmplot(data=df, color='black', s=4)
        else:
            fig = sns.swarmplot(data=df, x=xy[0], y=xy[1], color='black', s=4)
    if 'violin' in lower_plot_type:
        if xy == (None,None):
            fig = sns.violinplot(data=df)
        else:
            fig = sns.violinplot(data=df, x=xy[0], y=xy[1])
    if 'box' in lower_plot_type:
        if xy == (None,None):
            fig = sns.boxplot(data=df)
        else:
            fig = sns.boxplot(data=df, x=xy[0], y=xy[1])

    fig.yaxis.set_label_text(ylabel)
    fig.set_title(title)
    if xticks:
        fig.xaxis.set_ticklabels(xticks)
        fig.xaxis.set_label_text('')
        for tick in fig.xaxis.get_ticklabels():
            tick.set_fontsize(12)

    if pvalue:
        if xy==(None,None):
            _,pvalue = stats.ttest_ind(a=df.iloc[:,0], b=df.iloc[:,1])
            compare_tags = df.columns
        else:
            _,pvalue = stats.ttest_ind(a=df[df[xy[0]] == compare_tags[0]][xy[1]], b=df[df[xy[0]] == compare_tags[1]][xy[1]])
        fig.text(s='p-value = {:.03g}, {} v {}'.format(pvalue,compare_tags[0],compare_tags[1]), x=0, y=-.12, transform=fig.axes.transAxes, fontsize=12)
        
    sns.despine()
    plt.tight_layout()
    plt.savefig('{}{}.svg'.format(out,title.replace(' ','_')))
    plt.subplots_adjust(bottom=0.17, top=0.9)
    plt.savefig('{}{}.png'.format(out,title.replace(' ','_')), dpi=300)

    print('{}.png found in {}/'.format(title.replace(' ','_'), out))

def scatter_regression(df, s=150, alpha=0.3, line_color='dimgrey', svg=False, reg_stats=True, point_color='steelblue', title=None, 
                       xlabel=None,ylabel=None,IndexA=None,IndexB=None, annotate=None, Alabel='Group A', Blabel='Group B'):
    '''
    Scatter plot and Regression based on two matched vectors.
    Plots r-square and pvalue on .png

    Inputs
    ------
    df: dataframe to plot (column1 = x axis, column2= y axis)
    
    kwargs (defaults):
    s: point size (150)
    alpha: (0.3)
    line_color: regression line color (dimgrey)
    svg: make svg (False)
    stats: print R2 and pvalue on plot (True)
    point_color: (steelblue)
    title: string
    xlabel: string
    ylabel: string
    IndexA: set or list of genes to highlight red
    Alabel: string for IndexA group ('Group A')
    IndexB: set or list of genes to highlight blue
    annotate:  list of genes to annotate on the graph

    Returns
    -------
    None
    Prints file name and location
    Saves .png plot in scatter_regression/ folder in cwd with dpi=300.

    '''
    sns.set(context='paper',style="white", font_scale=3, font='Arial', 
            rc={"lines.linewidth": 2, 
                'figure.figsize': (9,9),  
                'font.size': 18, 'figure.dpi':300})
    fig,ax = plt.subplots()

    cols = df.columns.tolist()
    regplot = sns.regplot(x=cols[0], y=cols[1],data=df,scatter=True,
                          fit_reg = True, color=line_color,
                          scatter_kws= {'s': s, 'color': point_color, 'alpha':alpha}
                         )
    
    if xlabel:
        plt.xlabel(xlabel, labelpad = 10)
    if ylabel:
        plt.ylabel(ylabel, labelpad = 10)
    if title:
        regplot.set_title(title)
    if type(IndexA) in [list,set]:
        A = set(IndexA)
        Abool = [True if x in IndexA else False for x in df.index.tolist()]
        regplot = ax.scatter(df[Abool].iloc[:,0], df[Abool].iloc[:,1], marker = 'o', alpha = (alpha + .4 if alpha < .6 else 1), color='red', s=s, label= Alabel)
    if type(IndexB) in [list,set]:
        B = set(IndexB)
        Bbool = [True if x in IndexB else False for x in df.index.tolist()]
        regplot = ax.scatter(df[Bbool].iloc[:,0], df[Bbool].iloc[:,1], marker = 'o', alpha = (alpha + .3 if alpha < .7 else 1), color='mediumblue', s=s, label= Blabel)
    if type(annotate) in [list,set]:
        anno_df = df[[True if x in annotate else False for x in df.index.tolist()]]
        offx,offy = (df.iloc[:,:2].max()-df.iloc[:,:2].min())*.1
        for index, (x,y) in anno_df.iterrows():
            ax.annotate(index, xy=(x,y), xytext=((x-offx,y+offy) if y>=x else (x+offx, y-offy)), arrowprops={'arrowstyle':'-', 'color': 'black'})
    if reg_stats:
        r,pvalue = stats.pearsonr(x=df.iloc[:,0], y=df.iloc[:,1])
        ax.text(0,0 , 'r = {:.03g}; p-value = {:.03g}'.format(r,pvalue), fontsize = 25, transform=ax.transAxes)
    
    sns.despine(offset= 5)
    fig.tight_layout()

    os.makedirs('scatter_regression/', exist_ok=True)
    
    if svg:
        plt.savefig('scatter_regression/{}.svg'.format(title.replace(' ','_')))
    plt.savefig('scatter_regression/{}.png'.format(title.replace(' ','_')), dpi=300)

    print('{}.png found in {}/scatter_regression/'.format(title.replace(' ','_'), os.getcwd()))


def signature_heatmap(vst, sig, name, cluster_columns=False):
    '''
    Generate heatmap of differentially expressed genes using 
    variance stablized transfrmed log2counts.
    
    Inputs 
    ------
    vst = gene name is the index
    sig = set or list of signature
    name = name of file
    cluster_columns = bool (default = False)

    Outputs
    ------
    .png and .svg file of heatmap

    Returns
    -------
    None


    '''
    sns.set(font='Arial',font_scale=2, style='white', context='paper')
    vst['gene_name']=vst.index
    
    CM = sns.clustermap(vst[vst.gene_name.apply(lambda x: x in sig)].drop('gene_name',axis=1),
                        z_score=0, method='complete', cmap='RdBu_r', 
                        yticklabels=False, col_cluster=cluster_columns)
    CM.fig.suptitle(name)
    CM.savefig('{}_Heatmap.png'.format(name.replace(' ','_')), dpi=300)
    CM.savefig('{}_Heatmap.svg'.format(name.replace(' ','_')))

def image_display(file):
    display(Image(file))

def ssh_job(command_list, job_name, job_folder, project='nimerlab', threads=1, q='general', mem=3000):
    '''
    Sends job to LSF pegasus.ccs.miami.edu

    Inputs
    ------

    command_list: list of commands with new lines separated by commas
    job_name: string of job name (also used for log file)
    job_folder: string of folder to save err out and script files
    q: pegasus q, ie. 'bigmem', 'general' (default), 'parrallel'
    mem: integer memory requirement, default=3000 (3GB RAM)
    project: string pegasus project name (default = nimerlab)
    threads: integer of number of threads. default = 1
    ssh: whether or not to ssh into pegasus. Default=True

    Returns
    -------
    Tuple(rand_id, job_folder, prejob_files)

    '''
    
    jobfolder = job_folder if job_folder.endswith('/') else '{}/'.format(job_folder)

    os.system('ssh pegasus mkdir -p {}'.format(job_folder))

    rand_id = str(random.randint(0, 100000))
    str_comd_list =  '\n'.join(command_list)
    cmd = '''#!/bin/bash

#BSUB -J JOB_{job_name}_ID_{random_number}
#BSUB -R "rusage[mem={mem}]"
#BSUB -R "span[ptile={threads}]"
#BSUB -o {job_folder}{job_name_o}_logs_{rand_id}.stdout.%J
#BSUB -e {job_folder}{job_name_e}_logs_{rand_id}.stderr.%J
#BSUB -W 120:00
#BSUB -n 1
#BSUB -q {q}
#BSUB -P {project}

{str_comd_list}'''.format(job_name = job_name.replace(' ','_'),
                          job_folder=job_folder,
                          job_name_o=job_name.replace(' ','_'),
                          job_name_e=job_name.replace(' ','_'),
                          str_comd_list=str_comd_list,
                          random_number=rand_id,
                          rand_id=rand_id,
                          q=q,
                          mem=mem,
                          project=project,
                          threads=threads
                         )
    
    with open('{}.sh'.format(job_name.replace(' ','_')), 'w') as file:
        file.write(cmd)

    prejob_files = os.popen('ssh pegasus ls {}'.format(job_folder)).read().split('\n')[:-1]
    os.system('scp {}.sh pegasus:{}'.format(job_name.replace(' ','_'), job_folder))
    os.system('ssh pegasus "cd {}; bsub < {}.sh"'.format(job_folder, job_name.replace(' ','_')))
    print('Submitting {} as ID_{} from folder {}: {:%Y-%m-%d %H:%M:%S}'.format(job_name,rand_id, job_folder,datetime.now()))

    return (rand_id, job_folder, prejob_files, job_name)

def ssh_check(ID, job_folder, prejob_files=None, wait=True, return_filetype=None, load=False, check_IO_logs=None, sleep=10, job_name=''):
    '''
    Checks for pegasus jobs sent by ssh_job and prints contents of the log file.  
    Optionally copies and/or loads the results file.
    
    Inputs
    ------
    Job ID: Job ID
    wait: wait for processes to finish before returning, default=True
    job_folder: job folder to probe for results, (only if return_filetype specified)
    return_filetype: return file type (ex. .png will search for all .png in job_folder and import it)default=None
    display: whether or not to display imported file
    pre_list: list of contents of job folder brefore execution.
    check_IO_logs: read output from .err .out logs
    sleep: seconds to sleep (default 10)
    job_name: pepends local ssh folder with job name if provided

    Returns
    ------
    None

    '''
    jobfolder = job_folder if job_folder.endswith('/') else '{}/'.format(job_folder)

    jobs_list = os.popen('ssh pegasus bhist -w').read()
    job = [j for j in re.findall('ID_(\d+)', jobs_list) if j == ID]
    if len(job) != 0:
        print('Job ID_{} is not complete: {:%Y-%m-%d %H:%M:%S}'.format(ID, datetime.now()))
    else:
        if os.popen('''ssh pegasus "if [ -f {}/*_logs_{}.stderr* ]; then echo 'True' ; fi"'''.format(job_folder, ID)).read() == 'True\n':
            print('Job ID_{} is finished'.format(ID))
        else:
            print('There was likely an error in submission of Job ID_{}'.format(ID))

    if wait:
        running = True
        while running:
            jobs_list = os.popen('ssh pegasus "bhist -w"').read()
            job = [j for j in re.findall('ID_(\d+)', jobs_list) if j == ID]
            if len(job) == 0:
                running = False
            else:
                print('Waiting for jobs to finish... {:%Y-%m-%d %H:%M:%S}'.format(datetime.now()))
                time.sleep(sleep)
        print('Job ID_{} is finished'.format(ID))

    if load:
        os.makedirs('ssh_files/{}/'.format(ID), exist_ok=True)
        post_files = os.popen('ssh pegasus ls {}*{}'.format(job_folder,return_filetype)).read().split('\n')[:-1]

        if prejob_files is None:
            prejob_files = []
        import_files = [file for file in post_files if file not in prejob_files]

        for file in import_files:
            print('Copying {} to {}/ssh_files/{}{}/'.format(file, os.getcwd(), job_name,ID))
            os.system('scp pegasus:{} ssh_files/{}{}/{}'.format(file, job_name,ID,file.split('/')[-1]))
            image_display('ssh_files/{}{}/{}'.format(job_name,ID,file.split('/')[-1]))

    if check_IO_logs:
        logs= {'ErrorFile':'{}/*_logs_{}.stderr*'.format(job_folder, ID),
               'OutFile':'{}/*_logs_{}.stdout*'.format(job_folder, ID)
              }
        os.makedirs('logs/', exist_ok=True)
        for key,log in logs.items():
            os.system("scp 'pegasus:{}' 'logs/ID_{}_{}.txt'".format(log, ID, key))
            if os.path.isfile('logs/ID_{}_{}.txt'.format(ID,key)):
                print('logs/ID_{} {}:'.format(ID,key))
                with open('logs/ID_{}_{}.txt'.format(ID,key)) as file:
                    print(file.read())

def deeptools(regions, signals, matrix_name, out_name, pegasus_folder, title='', bps=(1500,1500,4000), type='center', scaled_names=('TSS','TES'), make=('matrix','heatmap','heatmap_group','profile', 'profile_group')):
    '''
    Inputs
    ------
    regions: dictionary {'region_name':'/path/to/ssh/bedfile'}
    signals: dictionary {'signal_name':'/path/to/ssh/bigwigfile'}
    matrix_name: string of matrix name or matrix to be named (before .matrix.gz)
    out_name: name for output file
    tite: plot title (optional)
    bps: tuple of region width on either side of center or scaled.  center ignores last number.  default is (1500,1500,4000)
    type: 'center' or 'scaled'
    scaled_names: optional names for scaled start and end (default ('TSS','TES'))
    make: tuple of deeptool commands.  options: matrix, heatmap, heatmap_group, profile, profile_group
    copy: bool.  Copy region and signal files to peagasus
    copy_folder: folder to copy into

    Returns
    -------
    string of commands for ssh_job

    '''
    pegasus_folder = pegasus_folder if pegasus_folder.endswith('/') else '{}/'.format(pegasus_folder)
    os.system("ssh pegasus 'mkdir {}'".format(pegasus_folder))

    make_lower = [x.lower() for x in make]

    if type.lower() == 'center':
        deepMat = 'reference-point --referencePoint center'
        deepHeat = "--refPointLabel 'Peak Center'"
        deepProf = "--refPointLabel 'Peak Center'"
    else:
        deepMat = 'scale-regions --regionBodyLength {bp3}'.format(str(bps[2]))
        deepHeat = '--startLabel {} --endLabel {}'.format(scaled_names[0],scaled_names[1])
        deepProf = '--startLabel {} --endLabel {}'.format(scaled_names[0],scaled_names[1])
        
    cmd_list = ['module rm python share-rpms65', 'source activate deeptools']
    
    print('Copying region files to pegasus...')
    for region in regions.values():
        if os.popen('''ssh pegasus "if [ -f {}{} ]; then echo 'True' ; fi"'''.format(pegasus_folder, region.split('/')[-1])).read() != 'True\n':
            print('Copying {} to pegasus at {}.'.format(region,pegasus_folder))
            os.system("scp {} pegasus:{}".format(region,pegasus_folder))
        else:
            print('{} found in {}.'.format(region,pegasus_folder))
    
    print('Copying signal files to pegasus...')
    for signal in signals.values():
        if os.popen('''ssh pegasus "if [ -f {}/{} ]; then echo 'True' ; fi"'''.format(pegasus_folder, signal.split('/')[-1])).read() != 'True\n':
            print('Copying {} to {}.'.format(signal, pegasus_folder))
            os.system("scp {} pegasus:{}".format(signal,pegasus_folder))

    pegasus_region_path = ' '.join(["{}{}".format(pegasus_folder,region_path.split('/')[-1]) for region_path in regions.values()])
    pegasus_signal_path = ' '.join(["{}{}".format(pegasus_folder,signal_path.split('/')[-1]) for signal_path in signals.values()])

    if 'matrix' in make_lower:
        computeMatrix = "computeMatrix {deepMat} -a {bp1} -b {bp2} -p 4 -R {region_path} -S {signal_path} --samplesLabel {signal_name} -o {matrix_name}.matrix.gz".format(
                                deepMat=deepMat,
                                bp1=str(bps[0]),
                                bp2=str(bps[1]),
                                region_path=pegasus_region_path,
                                signal_path=pegasus_signal_path,
                                signal_name=' '.join(["{}".format(signal_name) for signal_name in signals.keys()]),
                                matrix_name=matrix_name
                                )
        cmd_list.append(computeMatrix)

    if 'heatmap' in make_lower or 'heatmap_group' in make_lower:
        plotHeatmap_base = "plotHeatmap -m {matrix_name}.matrix.gz --dpi 300 {deepHeat} --regionsLabel {region_name} --plotTitle '{title}' --whatToShow 'heatmap and colorbar' --colorMap Reds -out {out_name}_heatmap".format(
                                matrix_name=matrix_name,
                                deepHeat=deepHeat,
                                region_name=' '.join(["{}".format(region_name) for region_name in regions.keys()]),
                                title=title,
                                out_name=out_name
                                )
        if 'heatmap' in make_lower:
            cmd_list.append("{}.png".format(plotHeatmap_base))
        if 'heatmap_group' in make_lower:
            cmd_list.append("{}_perGroup.png --perGroup".format(plotHeatmap_base))

    if 'profile' in make_lower or 'profile_group' in make_lower:
        plotProfile_base = "plotProfile -m {matrix_name}.matrix.gz --dpi 300 {deepProf} --plotTitle '{title}' --regionsLabel {region_name} -out {out_name}_profile".format(
                                matrix_name=matrix_name,
                                deepProf=deepProf,
                                title=title,
                                region_name = ' '.join(["{}".format(region_name) for region_name in regions.keys()]),
                                out_name=out_name
                                )
        if 'profile' in make_lower:
            cmd_list.append("{}.png".format(plotProfile_base))
        if 'profile_group' in make_lower:
            cmd_list.append("{}_perGroup.png --perGroup".format(plotProfile_base))

    return cmd_list


def gsea_barplot(out_dir,pos_file,neg_file,gmt_name,max_number=20):
    '''
    Inputs
    ------
    out_dir: directory output or '' for current directory
    pos_file: GSEA positive enrichment .xls file
    neg_file: GSEA negative enrichment .xls file
    gmt_name: name of enrichment (ex: Hallmarks)
    max_number: max number of significant sets to report (default 20)

    Returns
    -------
    string of save file

    '''
    
    out_dir = out_dir if out_dir.endswith('/') else '{}/'.format(out_dir)
    out_dir = '' if out_dir == '/' else out_dir 
    pos = pd.read_table(pos_file).head(max_number) if os.path.isfile(pos_file) else pd.DataFrame(columns=['FDR q-val'])
    pos[gmt_name] = [' '.join(name.split('_')[1:]) for name in pos.NAME.tolist()]
    neg = pd.read_table(neg_file).head(max_number) if os.path.isfile(neg_file) else pd.DataFrame(columns=['FDR q-val'])
    neg[gmt_name] = [' '.join(name.split('_')[1:]) for name in neg.NAME.tolist()]
    
    sns.set(context='paper', font='Arial',font_scale=.95, style='white', rc={'figure.dpi': 300, 'figure.figsize':(6,6)})
    fig,(ax1,ax2) = dk.plt.subplots(ncols=1, nrows=2)
    fig.suptitle('{} GSEA enrichment\n(q<0.05, max {})'.format(gmt_name, max_number))
    
    if len(pos[pos['FDR q-val'] < 0.05]) > 0:
        UP = sns.barplot(data=pos[pos['FDR q-val'] < 0.05], x = 'NES', y=gmt_name, color='firebrick', ax=ax1)
        UP.set_title('Positive Enrichment')
        sns.despine()
       
    if len(neg[neg['FDR q-val'] < 0.05]) > 0:
        DN = sns.barplot(data=neg[neg['FDR q-val'] < 0.05], x = 'NES', y=gmt_name, color='steelblue', ax=ax2)
        DN.set_title('Negative Enrichment')
        sns.despine()
    
    plt.tight_layout(h_pad=2,w_pad=1)
    plt.subplots_adjust(top=0.88)
    file='{}{}_GSEA_NES_plot.png'.format(out_dir,gmt_name)
    fig.savefig(file, dpi=300)
    plt.close()

    return file

'''
def get_text_positions(x_data, y_data, txt_width, txt_height):
    import numpy as np
    import matplotlib.pyplot as plt
    from numpy.random import *

    a = zip(y_data, x_data)
    text_positions = y_data.copy()
    for index, (y, x) in enumerate(a):
        local_text_positions = [i for i in a if i[0] > (y - txt_height) 
                            and (abs(i[1] - x) < txt_width * 2) and i != (y,x)]
        if local_text_positions:
            sorted_ltp = sorted(local_text_positions)
            if abs(sorted_ltp[0][0] - y) < txt_height: #True == collision
                differ = np.diff(sorted_ltp, axis=0)
                a[index] = (sorted_ltp[-1][0] + txt_height, a[index][1])
                text_positions[index] = sorted_ltp[-1][0] + txt_height
                for k, (j, m) in enumerate(differ):
                    #j is the vertical distance between words
                    if j > txt_height * 2: #if True then room to fit a word in
                        a[index] = (sorted_ltp[k][0] + txt_height, a[index][1])
                        text_positions[index] = sorted_ltp[k][0] + txt_height
                        break
    return text_positions

def text_plotter(x_data, y_data, text_positions, axis,txt_width,txt_height):
    import numpy as np
    import matplotlib.pyplot as plt
    from numpy.random import *
    
    for x,y,t in zip(x_data, y_data, text_positions):
        axis.text(x - txt_width, 1.01*t, '%d'%int(y),rotation=0, color='blue')
        if y != t:
            axis.arrow(x, t,0,y-t, color='red',alpha=0.3, width=txt_width*0.1, 
                       head_width=txt_width, head_length=txt_height*0.5, 
                       zorder=0,length_includes_head=True)


 
def TSS_enhancer_bar(tss_bed, enhancer_bed, bed_dict, title, folder):    
    ##

    import pandas as pd

    Numbers = dict(TSS=[], Enhancer=[], Other=[], Total=[], TSS_F=[], Enh_F=[], Other_F=[], name=[])
                   
    for name,bed in bed_dict.items():
        tss=len(bed.sort().merge().intersect(tss_bed.sort().merge()).intersect(enhancer_bed.sort().merge(),v=True))
        enhancer=len(bed.sort().merge().intersect(enhancer_bed.sort().merge()).intersect(tss_bed.sort().merge(),v=True))
        other = len(bed.sort().merge().intersect(enhancer_bed.sort().merge(), v=True).intersect(tss_bed.sort().merge(),v=True))
        total = tss + enhancer + other
        
        Numbers['TSS'].append(tss)
        Numbers['Enhancer'].append(enhancer)
        Numbers['Other'].append(other)
        Numbers['Total'].append(total)
        Numbers['TSS_F'].append(tss/total)
        Numbers['Enh_F'].append(enhancer/total)
        Numbers['Other_F'].append(other/total)
        Numbers['name'].append(name)
        
    df= pd.DataFrame(Numbers, index=Numbers['name'])
    df['Peak_type']=df.name.apply(lambda x: x.split(' ')[-1])
    df['IP']=df.name.apply(lambda x: x.split(' ')[0])
    df['Total_F']=1

    return df

def survival():
    ## import from lifelines

    ax=plt.subplot(111)

    kmf.fit(KM_vsd_top.Time, event_observed=KM_vsd_top.Event, label=['CARM1 High'])
    kmf.plot(show_censors=True,ax=ax)

    kmf.fit(KM_vsd_bottom.Time, event_observed=KM_vsd_bottom.Event, label=['CARM1 Lo'])
    kmf.plot(show_censors=True, ax=ax)

    plt.title("Survival using VSD Normalization")

    plt.savefig(Folder + "CARM1-vsd-Survival-CI.png", dpi=200)
    plt.savefig(Folder + "CARM1-vsd-Survival-CI.svg", dpi=200)

def surface_plot(df,name):
    #### import plotly

    lighting_effects = dict(ambient=0.4, diffuse=.9, roughness = 1, fresnel=5)
    data = [go.Surface(z=df.values.tolist(), colorscale='Reds', cmin=0, cmax=20, reversescale=False, lighting=lighting_effects)
           ]
    layout = go.Layout(
        width=800,
        height=700,
        autosize=False,
        title='{} Excess over Bliss'.format(name),
        scene=dict(
            xaxis=dict(
                title='PRMT5i',
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=False,
                ticktext=labels,
                tickvals=[i * 1 for i in range(len(labels))]
            ),
            yaxis=dict(
                title='PARPi',
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=False,
                ticktext=labels,
                tickvals=[i * 1 for i in range(len(labels))]
            ),
            zaxis=dict(
                title='eob',
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=False,
            ),
            aspectratio = dict( x=1, y=.7, z=.5 ),
            aspectmode = 'manual'
        )
    )

    fig = dict(data=data, layout=layout)
    plot(fig, filename='{}_eob.html'.format(name))

def chouchou(key, df):
    ## add options for changing names and inhibitors

    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt

    sns.set(context='paper', font_scale=2, font='Arial', style='white')
    plt.figure(dpi=300, figsize=(5,5))
    
    plt.plot([1,0], c='grey', linewidth=3)
    plt.scatter(x=df['PRMT5i(D/Dx)'], 
                y=df['PARPi(D/Dx)'], 
                c='red', s=120, 
                edgecolors='darkred', 
                linewidths=.7)
    plt.xlim(0,1.5)
    plt.ylim(0,1.5)
    plt.xticks(list(round(x*.2,1) for x in range(0,8)))
    plt.yticks(list(round(x*.2,1) for x in range(0,8)))
    plt.xlabel("D$_{PRMT5i}$/Dx$_{PRMT5i}$")
    plt.ylabel("D$_{PARPi}$/Dx$_{PARPi}$")
    plt.text(s=key,
             x=1.4, y=1.4, 
             horizontalalignment='right', 
             verticalalignment='top')

    plt.tight_layout()
    sns.despine()
    plt.savefig('{}_chou-chou_synergy_plot.png'.format(key), dpi=300)
    plt.savefig('{}_chou-chou_synergy_plot.svg'.format(key))
'''