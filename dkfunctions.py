#!/usr/bin/env python3

# various functions and mixins for downstream genomic and epigenomic anlyses

import os
import glob
import re
import random
from datetime import datetime
import time

from pybedtools import BedTool
import pandas as pd
import numpy as np

from tqdm import tqdm_notebook, tqdm

# Get Current Git Commit Hash for version
path = [x.replace(' ', r'\ ') for x in os.popen('echo $PYTHONPATH').read().split(':') if 'dkfunctions' in x.split('/')]

if len(path) > 0:
    version = os.popen(f'cd {path[0]}; git rev-parse HEAD').read()[:-1]
    __version__ = f'v0.1, Git SHA1: {version}'
else:
    __version__ = f'v0.1, {datetime.now():%Y-%m-%d}'


def val_folder(folder):
    folder = folder if folder.endswith('/') else f'{folder}/'
    folder = f'{os.getcwd()}/' if folder == '/' else folder
    os.makedirs(folder, exist_ok=True)
    return folder


def image_display(file):
    from IPython.display import Image, display
    display(Image(file))


def rplot(plot_func, filename, filetype, *args, **kwargs):
    from rpy2.robjects.packages import importr
    grdevices = importr('grDevices')
    filetype = filetype.lower()

    plot_types = {'png': grdevices.png,
                  'svg': grdevices.svg,
                  'pdf': grdevices.pdf
                  }

    plot_types[filetype](f'{filename}.{filetype}')
    return_object = plot_func(*args, **kwargs)
    grdevices.dev_off()

    if filetype == 'png':
        image_display(f'{filename}.{filetype}')

    return return_object


def read_pd(file, *args, **kwargs):
    if (file.split('.')[-1] == 'txt') or (file.split('.')[-1] == 'tab'):
        return pd.read_table(file, header=0, index_col=0, *args, **kwargs)
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
    print(x, file=open(f'{os.getcwd()}/R_out_{datetime.now():%Y-%m-%d}.txt', 'a'))


def alert_me(text):
    '''
    Send me a pop up alert to macosx.
    '''

    os.system(f'''osascript -e 'tell Application "System Events" to display dialog "{text}"' ''')


def tq_type():
    jupyter = True if os.environ['_'].endswith('jupyter') else False
    return tqdm_notebook if jupyter else tqdm


def peak_overlap_MC(df_dict, background, permutations=1000, seed=42, notebook=True):
    '''
    Monte Carlo simulation of peak overlaps in a given background
    pvalue calucated as liklihood over emperical random background overlap of shuffled peaks per chromosome.
    Inputs
    ------
    df_dict:  dictinoary of dataframes in bed format
    background genome space:  pybedtool bed of background genome space
    permutations:  number of permutations
    seed: random seed
    Returns
    -------
    pvalue
    '''

    np.random.seed(seed)
    tq = tq_type()

    # generate probability of chosing a chromosome region based on its size
    bregions = background.to_dataframe()
    bregions.index = range(len(bregions))
    bregions['Size'] = bregions.iloc[:, 2] - bregions.iloc[:, 1]
    total_size = bregions.Size.sum()
    bregions['fraction'] = bregions.Size / total_size

    bed_dict = {name: df.copy() for name, df in df_dict.items()}

    # determine length of each peak region
    for df in bed_dict.values():
        df['Length'] = df.iloc[:, 2] - df.iloc[:, 1]

    # determine baseline overlap intersect count of preshuffled peaks.
    A, B = bed_dict.values()

    overlap = len(BedTool.from_dataframe(A).sort().merge() + BedTool.from_dataframe(B).sort().merge())

    results = []

    for permutation in tq(range(permutations)):
        for df in bed_dict.values():
            # randomly pick a region in the background based on size distribution of the regions
            index_list = bregions.index.tolist()
            df_size = len(df)
            bregions_fraction = bregions.fraction

            first_pick = np.random.choice(index_list, size=df_size, p=bregions_fraction)
            lengths = df.Length.tolist()
            alternatives = np.random.choice(index_list, size=df_size, p=bregions_fraction)

            # repick regions if the peak length is larger than the region size (this part can be optimized)
            regions = []
            new_pick = 0
            for reg, length in zip(first_pick, lengths):
                reg_length = bregions.iloc[reg, 2] - bregions.iloc[reg, 1]
                if reg_length > length:
                    regions.append(reg)
                else:
                    while reg_length <= length:
                        new_reg = alternatives[new_pick]
                        reg_length = bregions.iloc[new_reg, 2] - bregions.iloc[new_reg, 1]
                        new_pick += 1
                    regions.append(new_reg)

            # assign the chromosome
            df.iloc[:, 0] = [bregions.iloc[x, 0] for x in regions]

            # randomly pick a start within the selected background region within the peak size constraints
            df.iloc[:, 1] = [np.random.randint(bregions.iloc[reg, 1], bregions.iloc[reg, 2] - length) for length, reg in zip(lengths, regions)]

            # assign end based on peak length
            df.iloc[:, 2] = df.iloc[:, 1] + df.Length

        new_overlap = len(BedTool.from_dataframe(A).sort().merge() + BedTool.from_dataframe(B).sort().merge())

        results.append(1 if new_overlap >= overlap else 0)

    p = (sum(results) + 1) / (len(results) + 1)

    A_name, B_name = df_dict.keys()

    print(f'Number of intersected peaks of {A_name} and {B_name}:  {overlap}')
    print(f'Number of times simulated intersections exceeded or equaled the actual overlap: {sum(results)}')
    print(f'Monte Carlo p-value estimate: {p}')

    return p


def gsea_dotplot(df_dict, title='', qthresh=0.05, top_term=None, gene_sets=[], dotsize_factor=4, figsize=(4, 10), out_dir='.'):
    '''
    Makes a dotplot of GSEA results with the dot size as the percent of genes in the leading edge and the color the NES.
    Plots only significant dots at given fdr theshold

    Inputs
    ------
    df_dict: dictionary of named GSEA results for the analysis. pandas df of gsea_report.xls (use pd.concat to combine pos and neg enrichments)
    name:  name used for title and filename
    qthresh: qvalue theshold for includsion
    pgene_sets:  list of gene sets to plot.  If empty, will plot all with FDR q value < 0.05
    top_term:  integer specifing top number of sets to plot (by qvalue).  None plots all.
    dot_size_factor:  scale to increase dot size for leading edge %
    out_dir: output directory

    Returns
    -------
    Gene_Sets used for plotting

    '''
    import matplotlib.pyplot as plt
    import seaborn as sns

    out_dir = val_folder(out_dir)

    index = []

    # get leading edge percentages
    for df in df_dict.values():
        if 'NAME' in df.columns.tolist():
            df.index = df.NAME
        df['le_tags'] = df['LEADING EDGE'].apply(lambda x: x.split('%')[0].split('=')[-1])
        df.sort_values(by='NES', ascending=False, inplace=True)
        index += df[df['FDR q-val'] < 0.05].index.tolist()

    index = list(set(index))

    # use gene_sets if provided
    if len(gene_sets) > 0:
        index = gene_sets

    # make master df
    data_df = pd.DataFrame()
    for name, df in df_dict.items():
        df['sample_name'] = name
        data_df = pd.concat([data_df, df.loc[index]])

    # extra filters
    data_df = data_df[data_df.sample_name.notna()]
    if top_term:
        index = list(set(data_df.sort_values(by='FDR q-val').head(top_term).index.tolist()))

    # reindex
    data_df['GS_NAME'] = data_df.index
    data_df.index = range(len(data_df))

    # make x coordinate
    samples = data_df.sample_name.unique()
    sample_number = len(samples)
    sample_x = {name: (x + .5) for name, x in zip(samples, range(sample_number))}
    data_df['x'] = data_df.sample_name.map(sample_x)

    # make y coordinate
    gene_set = list(index[::-1])
    gene_set_number = len(gene_set)
    sample_y = {name: y for name, y in zip(gene_set, range(gene_set_number))}
    data_df['y'] = data_df.GS_NAME.map(sample_y)

    # filter for significance and make dot size from leading edge percentage
    data_df['sig_tags'] = data_df[['FDR q-val', 'le_tags']].apply(lambda x: 0 if float(x[0]) > qthresh else float(x[1]), axis=1)
    data_df['area'] = data_df['sig_tags'] * dotsize_factor

    plot_df = data_df[data_df.GS_NAME.isin(index)].copy()

    # plot
    plt.clf()
    sns.set(context='paper', style='white', font='Arial', rc={'figure.dpi': 300})

    fig, ax = plt.subplots(figsize=figsize)
    sc = ax.scatter(x=plot_df.x, y=plot_df.y, s=plot_df.area, edgecolors='face', c=plot_df.NES, cmap='RdBu_r')

    # format y axis
    ax.yaxis.set_major_locator(plt.FixedLocator(plot_df.y))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(plot_df.GS_NAME))
    ax.set_yticklabels(plot_df.GS_NAME.apply(lambda x: x.replace('_', ' ')), fontsize=16)

    # format x axis
    ax.set_xlim(0, sample_number)
    ax.xaxis.set_major_locator(plt.FixedLocator(plot_df.x))
    ax.xaxis.set_major_formatter(plt.FixedFormatter(plot_df.sample_name))
    ax.set_xticklabels(plot_df.sample_name, fontsize=16, rotation=45)

    # add colorbar
    cax = fig.add_axes([0.95, 0.20, 0.03, 0.22])
    cbar = fig.colorbar(sc, cax=cax,)
    cbar.ax.tick_params(right=True)
    cbar.ax.set_title('NES', loc='left', fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    # add legend
    markers = []
    min_value = plot_df[plot_df.sig_tags > 0].sig_tags.min()
    max_value = plot_df.sig_tags.max()

    rounded_min = int(10 * round((min_value - 5) / 10))
    rounded_max = int(10 * round((max_value + 5) / 10))  # rounds up to nearest ten (ie 61 --> 70)

    sizes = [x for x in range(rounded_min, rounded_max + 1, 10)]
    for size in sizes:
        markers.append(ax.scatter([], [], s=size * dotsize_factor, c='k'))
    legend = ax.legend(markers, sizes, prop={'size': 12})
    legend.set_title('Leading Edge (%)', prop={'size': 12})

    # offset legend
    bb = legend.get_bbox_to_anchor().inverse_transformed(ax.transAxes)
    xOffset = .6
    yOffset = 0
    bb.x0 += xOffset
    bb.x1 += xOffset
    bb.y0 += yOffset
    bb.y1 += yOffset
    legend.set_bbox_to_anchor(bb, transform=ax.transAxes)

    # set title
    ax.set_title(title.replace('_', ' '), fontsize=20)

    sns.despine()

    fig.savefig(f'{out_dir}{title.replace(" ", "_")}.png', bbox_inches='tight')
    fig.savefig(f'{out_dir}{title.replace(" ", "_")}.svg', bbox_inches='tight')

    return plot_df


def annotate_peaks(dict_of_dfs, folder, genome, db='UCSC', check=False, TSS=[-3000,3000], clean=False):
    '''
    Annotate a dictionary of dataframes from bed files to the genome using ChIPseeker and Ensembl annotations.
    Inputs
    ------
    dict_of_beds: dictionary of bed files
    folder: output folder
    genome: hg38, hg19, mm10
    db: default UCSC, but can also accept Ensembl
    TSS: list of regions around TSS to annotate as promoter
    check: bool. checks whether annotation file already exists
    Returns
    -------
    dictionary of annotated bed files as dataframe
    '''
    import rpy2.robjects as ro
    import rpy2.rinterface as ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri

    pandas2ri.activate()

    tq = tq_type()

    ri.set_writeconsole_regular(rout_write)
    ri.set_writeconsole_warnerror(rout_write)

    chipseeker = importr('ChIPseeker')
    genomicFeatures = importr('GenomicFeatures')
    makeGR = ro.r("makeGRangesFromDataFrame")
    as_df = ro.r("as.data.frame")

    check_df = {key: os.path.isfile(f'{folder}{key.replace(" ","_")}_annotated.txt') for key in dict_of_dfs.keys()}
    return_bool = False not in set(check_df.values())
    if return_bool & check:
        return {f'{key}_annotated': pd.from_csv(f'{folder}{key.replace(" ","_")}_annotated.txt', index_col=0, header=0, sep="\t") for key in dict_of_dfs.keys()}

    species = ('Mmusculus' if genome.lower() == 'mm10' else 'Hsapiens')
    if db.lower() == 'ucsc':
        TxDb = importr(f'TxDb.{species}.UCSC.{genome.lower()}.knownGene')
        txdb = ro.r(f'txdb <- TxDb.{species}.UCSC.{genome.lower()}.knownGene')
    elif db.lower() == 'ensembl':
        TxDb = importr(f'TxDb.{species}.UCSC.{genome.lower()}.ensGene')
        txdb = ro.r(f'txdb <- TxDb.{species}.UCSC.{genome.lower()}.ensGene')
    else:
        raise ValueError('UCSC or Ensembl only.')

    os.makedirs(folder, exist_ok=True)

    if genome.lower() == 'mm10':
        annoDb = importr('org.Mm.eg.db')
        anno = 'org.Mm.eg.db'
    elif genome.lower() == 'hg38' or genome.lower() == 'hg19':
        annoDb = importr('org.Hs.eg.db')
        anno = 'org.Hs.eg.db'

    return_dict = {}

    print('Annotating Peaks...')
    for key, df in tq(dict_of_dfs.items()):
        if check & check_df[key]:
            return_dict[f'{key}_annotated'] = pd.from_csv(f'{folder}{key.replace(" ","_")}_annotated.txt', index_col=0, header=0, sep="\t")
        else:
            col_len = len(df.columns)
            df.columns = ["chr", "start", "end"] + list(range(col_len - 3))
            GR = makeGR(df)
            GR_anno = chipseeker.annotatePeak(GR, overlap='TSS', TxDb=txdb, annoDb=anno, tssRegion=ro.IntVector(TSS)) #switched to TSS on 10/02/2019
            return_dict[f'{key}_annotated'] = ro.pandas2ri.ri2py(chipseeker.as_data_frame_csAnno(GR_anno))
            return_dict[f'{key}_annotated'].to_excel(f'{folder}{key.replace(" ","_")}_annotated.xlsx')

    if clean:
        for k,df in return_dict.items():
            df['Anno'] = df.annotation.apply(lambda x: 'Promoter' if x.split(' ')[0] == 'Promoter' else x)
            df['Anno'] = df.Anno.apply(lambda x: 'Intergenic' if x.split(' ')[0] in ['Downstream', 'Distal'] else x)
            df['Anno'] = df.Anno.apply(lambda x: x.split(' ')[0] if x.split(' ')[0] in ['Intron', 'Exon'] else x)

    return return_dict


"""
def plot_peak_genomic_annotation(dict_of_df, folder, genome):
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
    """


def plot_venn2(Series, string_name_of_overlap, folder):
    '''
    Series with with overlaps 10,01,11
    Plots a 2 way venn.
    Saves to file.
    '''
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib_venn import venn2, venn2_circles

    folder = f'{folder}venn2/' if folder.endswith('/') else f'{folder}/venn2/'
    os.makedirs(folder, exist_ok=True)

    plt.figure(figsize=(7, 7))

    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
            }

    plt.rc('font', **font)

    # make venn
    sns.set(style='white', font='Arial')
    venn_plot = venn2(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]), set_labels=[name.replace('_', ' ') for name in Series.index.tolist()])
    patch = ['10', '01', '11']
    colors = ['green', 'blue', 'teal']
    for patch, color in zip(patch, colors):
        venn_plot.get_patch_by_id(patch).set_color('none')
        venn_plot.get_patch_by_id(patch).set_alpha(.4)
        venn_plot.get_patch_by_id(patch).set_edgecolor('none')

    c = venn2_circles(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]))
    colors_test = ['green', 'blue']
    for circle, color in zip(c, colors_test):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(3)

    plt.title(string_name_of_overlap.replace('_', ' ') + " overlaps")
    plt.tight_layout()
    name = string_name_of_overlap.replace('_', ' ').replace('\n', '_')
    plt.savefig(f"{folder}{name}-overlap.svg")
    plt.savefig(f"{folder}{name}-overlap.png", dpi=300)
    plt.close()

    image_display(f"{folder}{name}-overlap.png")


def plot_venn2_set(dict_of_sets, string_name_of_overlap, folder, pvalue=False, total_genes=None):
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
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib_venn import venn2, venn2_circles
    from scipy import stats

    folder = f'{folder}venn2/' if folder.endswith('/') else f'{folder}/venn2/'
    os.makedirs(folder, exist_ok=True)

    plt.figure(figsize=(7, 7))

    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
            }

    plt.rc('font', **font)

    set_list = []
    set_names = []
    for name, setlist in dict_of_sets.items():
        set_list.append(setlist)
        set_names.append(name.replace('_', ' '))

    # make venn
    sns.set(style='white', font='Arial')
    venn_plot = venn2(subsets=set_list, set_labels=set_names)
    patch = ['10', '01', '11']
    colors = ['green', 'blue', 'teal']
    for patch, color in zip(patch, colors):
        venn_plot.get_patch_by_id(patch).set_color('none')
        venn_plot.get_patch_by_id(patch).set_alpha(.4)
        venn_plot.get_patch_by_id(patch).set_edgecolor('none')

    c = venn2_circles(subsets=set_list)
    colors_test = ['green', 'blue']
    for circle, color in zip(c, colors_test):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(3)

    if None not in [pvalue, total_genes]:
        intersection_N = len(set_list[0] & set_list[1])
        pvalue = stats.hypergeom.sf(intersection_N, total_genes, len(set_list[0]), len(set_list[1]))
        pvalue_string = f'= {pvalue:.03g}' if pvalue > 1e-5 else '< 1e-5'
        plt.text(0, -.05, f'p-value {pvalue_string}', fontsize=10, transform=c[1].axes.transAxes)

    plt.title(string_name_of_overlap.replace('_', ' ') + " overlaps")
    plt.tight_layout()
    plt.savefig(f"{folder}{string_name_of_overlap.replace(' ', '_')}-overlap.svg")
    plt.savefig(f"{folder}{string_name_of_overlap.replace(' ', '_')}-overlap.png", dpi=300)

    plt.close()
    image_display(f"{folder}{string_name_of_overlap.replace(' ', '_')}-overlap.png")


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
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib_venn import venn3, venn3_circles

    folder = f'{folder}venn3/' if folder.endswith('/') else f'{folder}/venn3/'
    os.makedirs(folder, exist_ok=True)

    plt.clf()
    sns.set(style='white', context='paper', font_scale=2, rc={'figure.figsize': (7, 7)})

#    font = {'family': 'sans-serif',
#            'weight': 'normal',
#            'size': 16,
#            }

#   plt.rc('font', **font)

    set_list = []
    set_names = []
    for name, setlist in dict_of_sets.items():
        set_list.append(setlist)
        set_names.append(name.replace('_', ' '))

    # make venn
    venn_plot = venn3(subsets=set_list, set_labels=set_names)
    patch = ['100', '110', '101', '010', '011', '001', '111']
    for p in patch:
        if venn_plot.get_patch_by_id(p):
            venn_plot.get_patch_by_id(p).set_color('none')
            venn_plot.get_patch_by_id(p).set_alpha(.4)
            venn_plot.get_patch_by_id(p).set_edgecolor('none')

    # make
    c = venn3_circles(subsets=set_list)
    colors_list = ['green', 'blue', 'grey']
    for circle, color in zip(c, colors_list):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(3)

    plt.title(f"{string_name_of_overlap.replace('_', ' ')} Overlaps")
    plt.tight_layout()
    plt.savefig(f"{folder}{string_name_of_overlap.replace(' ','_')}-overlap.svg")
    plt.savefig(f"{folder}{string_name_of_overlap.replace(' ','_')}-overlap.png", dpi=300)

    plt.close()
    image_display(f"{folder}{string_name_of_overlap.replace(' ','_')}-overlap.png")


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
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib_venn import venn3, venn3_circles

    folder = f'{folder}venn3/' if folder.endswith('/') else f'{folder}/venn3/'
    os.makedirs(folder, exist_ok=True)

    plt.clf()

    sns.set(style='white', context='paper', font_scale=1, rc={'figure.figsize': (7, 7)})

    # font = {'family': 'sans-serif',
    #        'weight': 'normal',
    #        'size': 16,
    #        }

    # plt.rc('font', **font)

    # make venn
    venn_plot = venn3(subsets=element_list, set_labels=[name.replace('_', ' ') for name in set_labels])
    patch = ['100', '110', '101', '010', '011', '001', '111']
    for p in patch:
        if venn_plot.get_patch_by_id(p):
            venn_plot.get_patch_by_id(p).set_color('none')
            venn_plot.get_patch_by_id(p).set_alpha(.4)
            venn_plot.get_patch_by_id(p).set_edgecolor('none')

    # make
    c = venn3_circles(subsets=element_list)
    colors_list = ['green', 'blue', 'grey']
    for circle, color in zip(c, colors_list):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(3)

    plt.title(f"{string_name_of_overlap.replace('_', ' ')} Overlaps")
    plt.tight_layout()
    plt.savefig(f"{folder}{string_name_of_overlap.replace(' ', '_')}-overlap.svg")
    plt.savefig(f"{folder}{string_name_of_overlap.replace(' ', '_')}-overlap.png", dpi=300)

    plt.close()
    image_display(f"{folder}{string_name_of_overlap.replace(' ', '_')}-overlap.png")


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

    Folder = f'{os.getcwd()}/'
    subfolder = f"{names[0].replace(' ', '_')}_{names[1].replace(' ', '_')}_overlap/"

    out = f'{Folder}{subfolder}'
    os.makedirs(out, exist_ok=True)
    print(f'Output files are found in {out}')

    masterfile = bed_dict[names[0]].cat(bed_dict[names[1]]).sort().merge()
    sorted_dict = {key: bed.sort().merge() for key, bed in bed_dict.items()}
    overlap_dict = {'overlap': masterfile.intersect(sorted_dict[names[0]]).intersect(sorted_dict[names[1]])}
    for key, bed in sorted_dict.items():
        other = {other_key: other_bed for other_key, other_bed in sorted_dict.items() if other_key != key}
        overlap_dict['{}_unique_peak'.format(key)] = masterfile.intersect(sorted_dict[key]).intersect(list(other.values())[0], v=True)

    for key, bed in overlap_dict.items():
        bed.to_dataframe().to_csv('{}{}{}-unique-peaks-from-mergedPeaks.bed'.format(Folder, subfolder, key.replace(' ', '_')),
                                  header=None, index=None, sep="\t")

    overlap_numbers = pd.Series({names[0]: len(overlap_dict['{}_unique_peak'.format(names[0])]),
                                 names[1]: len(overlap_dict['{}_unique_peak'.format(names[1])]),
                                 'overlap': len(overlap_dict['overlap'])
                                 },
                                index=[names[0], names[1], 'overlap']
                                )

    # Venn
    plot_venn2(overlap_numbers,
               '{} and\n{} peak'.format(names[0], names[1]),
               '{}{}'.format(Folder, subfolder)
               )

    if bool(genome):
        print('Annotating overlaping peaks...')
        # Annotate with ChIPseeker
        unikey = '{}_unique'
        unianno = '{}_unique_annotated'
        return_dict = annotate_peaks({unikey.format(key): bed.to_dataframe() for key, bed in overlap_dict.items()}, '{}{}'.format(Folder, subfolder), genome=genome)

        Set1_unique = set(return_dict[unianno.format('{}_unique_peak'.format(names[0]))].SYMBOL.unique().tolist())
        Set2_unique = set(return_dict[unianno.format('{}_unique_peak'.format(names[1]))].SYMBOL.unique().tolist())
        Overlap_Set = set(return_dict[unianno.format('overlap')].SYMBOL.unique().tolist())

        venn2_dict = {names[0]: (Set1_unique | Overlap_Set),
                      names[1]: (Set2_unique | Overlap_Set)
                      }

        plot_venn2_set(venn2_dict,
                       '{} and {}\nannotated gene'.format(names[0], names[1]),
                       '{}{}'.format(Folder, subfolder)
                       )

        gene_overlaps = {}
        gene_overlaps['{}_unique_genes'.format(names[0])] = Set1_unique - (Set2_unique | Overlap_Set)
        gene_overlaps['{}_unique_genes'.format(names[1])] = Set2_unique - (Set1_unique | Overlap_Set)
        gene_overlaps['Overlap_Gene_Set'] = (Set1_unique & Set2_unique) | Overlap_Set

        for key, gene_set in gene_overlaps.items():
            with open(f'{Folder}{subfolder}{key}.txt', 'w') as file:
                for gene in gene_set:
                    file.write(f'{gene}\n')

        for key, item in gene_overlaps.items():
            return_dict[key] = item

        for key, df in overlap_dict.items():
            return_dict[key] = df

    else:
        return_dict = overlap_dict

    return return_dict


def enrichr(gene_list, description, out_dir, scan=None, max_terms=10, load=False, figsize=(12, 6), dotplot=False, save_type='png'):
    '''
    Performs GO Molecular Function, GO Biological Process and KEGG enrichment on a gene list.
    Uses enrichr.
    Inputs
    ------
    gene_list: list of genes to perform enrichment on
    description: string description for title
    out_dir: output director
    scan: dictionary with additional enrichr dbs to scan (http://amp.pharm.mssm.edu/Enrichr/#stats)
    max_terms: limit return plot to this max
    load: load results
    figsize: change fig size
    Returns
    -------
    None
    '''
    import gseapy

    out_dir = val_folder(out_dir)
    tq = tq_type()

    testscan = {'KEGG': 'KEGG_2016',
                'GO_biological_process': 'GO_Biological_Process_2017b',
                'ChIP-X_Consensus_TFs': 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
                'ChEA': 'ChEA_2016',
                'OMIM_Disease': 'OMIM_Disease'
                }

    if isinstance(scan, dict):
        testscan = {**testscan, **scan}

    for nick, name in tq(testscan.items()):
        try:
            res = gseapy.enrichr(gene_list=gene_list,
                                 figsize=figsize,
                                 top_term=max_terms,
                                 description=f'{description}_{nick}',
                                 gene_sets=name,
                                 outdir=out_dir,
                                 format=save_type
                                 )

            if dotplot:
                dtplt = f'{out_dir}{description}_{name}_dotplot.{save_type}'
                gseapy.dotplot(res,
                               title=f'{description} {name} DotPlot',
                               top_term=20,
                               ofname=dtplt
                               )

            if load & (save_type == 'png'):
                png = f'{out_dir}{name}.{description}_{nick}.enrichr.reports.png'
                if os.path.isfile(png):
                    print(f'{description}_{nick}:')
                    image_display(png)
                if dotplot:
                    image_display(dtplt)
        except:
            print('Error with enrichr')

    with open(f'{out_dir}{description}_genes.txt', 'w') as fp:
        for gene in set(gene_list):
            fp.write(f'{gene}\n')


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

    Folder = f'{os.getcwd()}/'
    subfolder = f"{names[0].replace(' ', '_')}-{ names[1].replace(' ', '_')}-{names[2].replace(' ', '_')}-overlap/"

    out = f'{Folder}{subfolder}'
    os.makedirs(out, exist_ok=True)
    print(f'Output files are found in {out}')
    print(f'A: {names[0]}, B: {names[1]}, C: {names[2]}')

    master = bed_dict[names[0]].cat(bed_dict[names[1]]).cat(bed_dict[names[2]]).sort().merge()

    A = bed_dict[names[0]].sort().merge()
    B = bed_dict[names[1]].sort().merge()
    C = bed_dict[names[2]].sort().merge()

    sorted_dict = OrderedDict({'master': master, 'A': A, 'B': B, 'C': C})
    sorted_dict['A_bc'] = (master + A - B - C)
    sorted_dict['aB_c'] = (master + B - A - C)
    sorted_dict['A_B_c'] = (master + A + B - C)
    sorted_dict['abC_'] = (master + C - A - B)
    sorted_dict['A_bC_'] = (master + A + C - B)
    sorted_dict['aB_C_'] = (master + B + C - A)
    sorted_dict['A_B_C_'] = (master + A + B + C)

    labTup = tuple(key for key in sorted_dict.keys())
    lenTup = tuple(len(bed) for bed in sorted_dict.values())

    print(f'{labTup}\n{lenTup}')

    plot_venn3_counts(lenTup[4:], names, f"{'_'.join(names)}-peak-overlaps", out)

    for key, bed in sorted_dict.items():
        if len(bed) > 1:
            bed.to_dataframe().to_csv(f'{out}{key.replace(" ", "_")}-peaks-from-mergedPeaks.bed', header=None, index=None, sep="\t")

    if bool(genome):
        print('Annotating ovelapped peaks...')
        unikey = '{}'
        unianno = '{}_annotated'
        return_dict = annotate_peaks({unikey.format(key): bed.to_dataframe() for key, bed in sorted_dict.items()}, out, genome=genome)

        Set1 = set(return_dict[unianno.format('A')].SYMBOL.unique().tolist())
        Set2 = set(return_dict[unianno.format('B')].SYMBOL.unique().tolist())
        Set3 = set(return_dict[unianno.format('C')].SYMBOL.unique().tolist())

        plot_venn3_set({names[0]: Set1, names[1]: Set2, names[2]: Set3}, f'{names[0]}_{names[1]}_{names[2]}-gene-overlaps', out)

    return sorted_dict if genome is None else {**sorted_dict, **return_dict}


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
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set(context='paper', font='Arial', style='white', font_scale=2)

    plot = sns.barplot(x=x, y=y, data=data)
    plot.set_title(title.replace('_', ' '))
    plot.set_ylabel('')

    sns.despine()
    sns.utils.plt.savefig('{}.png'.format(title.replace(' ', '_')), dpi=300)

    plt.close()
    image_display('{}.png'.format(title.replace(' ', '_')))


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

    out_dir = '{pwd}/{name}/'.format(pwd=os.getcwd(), name=name.replace(' ', '_'))
    os.makedirs(out_dir, exist_ok=True)

    count = 0
    for key, genes in dict_of_sets.items():
        count = max(count, len(genes))

    df = pd.DataFrame(index=range(1, count + 1))

    for key, genes in dict_of_sets.items():
        df[key] = pd.Series(list(genes) + ['NA'] * (count - len(genes)))

    df.to_excel('{}/{}.xls'.format(out_dir, name.replace(' ', '_')), index=False)

    return df


def plot_col(df, title, ylabel, out='', xy=(None, None), xticks=[''], plot_type=['violin'], pvalue=False, compare_tags=None):
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
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats

    out = val_folder(out)

    plt.clf()
    sns.set(context='paper', font='Arial', font_scale=2, style='white', rc={'figure.dpi': 300, 'figure.figsize': (5, 6)})

    if type(plot_type) != list:
        plot_type = plot_type.split()
    lower_plot_type = [x.lower() for x in plot_type]

    if len(lower_plot_type) == 0:
        raise IOError('Input a plot type.')
    elif True not in {x in lower_plot_type for x in ['violin', 'box', 'swarm']}:
        raise IOError('Did not recognize plot type.')

    if 'swarm' in lower_plot_type:
        if xy == (None, None):
            fig = sns.swarmplot(data=df, color='black', s=4)
        else:
            fig = sns.swarmplot(data=df, x=xy[0], y=xy[1], color='black', s=4)
    if 'violin' in lower_plot_type:
        if xy == (None, None):
            fig = sns.violinplot(data=df)
        else:
            fig = sns.violinplot(data=df, x=xy[0], y=xy[1])
    if 'box' in lower_plot_type:
        if xy == (None, None):
            fig = sns.boxplot(data=df)
        else:
            fig = sns.boxplot(data=df, x=xy[0], y=xy[1])

    fig.yaxis.set_label_text(ylabel)
    fig.set_title(title.replace('_', ' '))
    if xticks:
        fig.xaxis.set_ticklabels(xticks)
        fig.xaxis.set_label_text('')
        for tick in fig.xaxis.get_ticklabels():
            tick.set_fontsize(12)

    if pvalue:
        if xy == (None, None):
            _, pvalue = stats.ttest_ind(a=df.iloc[:, 0], b=df.iloc[:, 1])
            compare_tags = df.columns
        else:
            _, pvalue = stats.ttest_ind(a=df[df[xy[0]] == compare_tags[0]][xy[1]], b=df[df[xy[0]] == compare_tags[1]][xy[1]])
        fig.text(s='p-value = {:.03g}, {} v {}'.format(pvalue, compare_tags[0], compare_tags[1]), x=0, y=-.12, transform=fig.axes.transAxes, fontsize=12)

    sns.despine()
    plt.tight_layout()
    plt.savefig('{}{}.svg'.format(out, title.replace(' ', '_')))
    plt.subplots_adjust(bottom=0.17, top=0.9)
    plt.savefig('{}{}.png'.format(out, title.replace(' ', '_')), dpi=300)

    print('{}.png found in {}/'.format(title.replace(' ', '_'), out))
    plt.close()
    image_display('{}{}.png'.format(out, title.replace(' ', '_')))


def scatter_regression(df, s=150, alpha=0.3, line_color='dimgrey', svg=False, reg_stats=True, point_color='steelblue', title=None,
                       xlabel=None, ylabel=None, IndexA=None, IndexB=None, annotate=None, Alabel='Group A', Blabel='Group B'):
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
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats

    sns.set(context='paper', style="white", font_scale=3, font='Arial',
            rc={"lines.linewidth": 2,
                'figure.figsize': (9, 9),
                'font.size': 18, 'figure.dpi': 300})
    fig, ax = plt.subplots()

    cols = df.columns.tolist()
    regplot = sns.regplot(x=cols[0], y=cols[1], data=df, scatter=True,
                          fit_reg=True, color=line_color,
                          scatter_kws={'s': s, 'color': point_color, 'alpha': alpha}
                          )

    if xlabel:
        plt.xlabel(xlabel, labelpad=10)
    if ylabel:
        plt.ylabel(ylabel, labelpad=10)
    if title:
        regplot.set_title(title.replace('_', ' '))
    if type(IndexA) in [list, set]:
        # A = set(IndexA)
        Abool = [True if x in IndexA else False for x in df.index.tolist()]
        regplot = ax.scatter(df[Abool].iloc[:, 0], df[Abool].iloc[:, 1], marker='o', alpha=(alpha + .4 if alpha < .6 else 1), color='red', s=s, label=Alabel)
    if type(IndexB) in [list, set]:
        # B = set(IndexB)
        Bbool = [True if x in IndexB else False for x in df.index.tolist()]
        regplot = ax.scatter(df[Bbool].iloc[:, 0], df[Bbool].iloc[:, 1], marker='o', alpha=(alpha + .3 if alpha < .7 else 1), color='mediumblue', s=s, label=Blabel)
    if type(annotate) in [list, set]:
        anno_df = df[[True if x in annotate else False for x in df.index.tolist()]]
        offx, offy = (df.iloc[:, :2].max() - df.iloc[:, :2].min()) * .1
        for index, (x, y) in anno_df.iterrows():
            ax.annotate(index, xy=(x, y), xytext=((x - offx, y + offy) if y >= x else (x + offx, y - offy)), arrowprops={'arrowstyle': '-', 'color': 'black'})
    if reg_stats:
        r, pvalue = stats.pearsonr(x=df.iloc[:, 0], y=df.iloc[:, 1])
        ax.text(0, 0, 'r = {:.03g}; p-value = {:.03g}'.format(r, pvalue), fontsize=25, transform=ax.transAxes)

    sns.despine(offset=5)
    fig.tight_layout()

    os.makedirs('scatter_regression/', exist_ok=True)

    if svg:
        plt.savefig('scatter_regression/{}.svg'.format(title.replace(' ', '_')))
    plt.savefig('scatter_regression/{}.png'.format(title.replace(' ', '_')), dpi=300)

    print('{}.png found in {}/scatter_regression/'.format(title.replace(' ', '_'), os.getcwd()))
    plt.close()
    image_display('scatter_regression/{}.png'.format(title.replace(' ', '_')))


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
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set(font='Arial', font_scale=2, style='white', context='paper')
    vst['gene_name'] = vst.index

    CM = sns.clustermap(vst[vst.gene_name.apply(lambda x: x in sig)].drop('gene_name', axis=1),
                        z_score=0, method='complete', cmap='RdBu_r',
                        yticklabels=False, col_cluster=cluster_columns)
    CM.fig.suptitle(name.replace('_', ' '))
    CM.savefig('{}_Heatmap.png'.format(name.replace(' ', '_')), dpi=300)
    CM.savefig('{}_Heatmap.svg'.format(name.replace(' ', '_')))

    plt.close()
    image_display('{}_Heatmap.png'.format(name.replace(' ', '_')))


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

    job_folder = job_folder if job_folder.endswith('/') else f'{job_folder}/'

    os.system(f'ssh pegasus mkdir -p {job_folder}')

    rand_id = str(random.randint(0, 100000))
    str_comd_list = '\n'.join(command_list)
    cmd = '\n'.join(['#!/bin/bash',
                     '',
                     f"#BSUB -J ID_{rand_id}_JOB_{job_name.replace(' ','_')}",
                     f'#BSUB -R "rusage[mem={mem}]"',
                     f'#BSUB -R "span[ptile={threads}]"',
                     f"#BSUB -o {job_folder}{job_name.replace(' ','_')}_logs_{rand_id}.stdout.%J",
                     f"#BSUB -e {job_folder}{job_name.replace(' ','_')}_logs_{rand_id}.stderr.%J",
                     '#BSUB -W 120:00',
                     f'#BSUB -n {threads}',
                     f'#BSUB -q {q}',
                     f'#BSUB -P {project}',
                     '',
                     f'{str_comd_list}'
                     ])

    with open(f'{job_name.replace(" ","_")}.sh', 'w') as file:
        file.write(cmd)

    prejob_files = os.popen(f'ssh pegasus ls {job_folder}').read().split('\n')[:-1]
    os.system(f'''ssh pegasus "mkdir -p {job_folder}"''')
    os.system(f'scp {job_name.replace(" ", "_")}.sh pegasus:{job_folder}')
    os.system(f'''ssh pegasus "cd {job_folder}; bsub < {job_name.replace(' ','_')}.sh"''')
    print(f'Submitting {job_name} as ID_{rand_id} from folder {job_folder}: {datetime.now():%Y-%m-%d %H:%M:%S}')

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
    job_folder = val_folder(job_folder)
    jobs_list = os.popen('ssh pegasus bhist -w').read()
    job = [j for j in re.findall(r'ID_(\d+)', jobs_list) if j == ID]
    if len(job) != 0:
        print(f'Job ID_{ID} is not complete: {datetime.now():%Y-%m-%d %H:%M:%S}')
    else:
        if os.popen('''ssh pegasus "if [ -f {}/*_logs_{}.stderr* ]; then echo 'True' ; fi"'''.format(job_folder, ID)).read() == 'True\n':
            print(f'Job ID_{ID} is finished')
        else:
            print(f'There was likely an error in submission of Job ID_{ID}')

    if wait:
        running = True
        while running:
            jobs_list = os.popen('ssh pegasus "bhist -w"').read()
            job = [j for j in re.findall(r'ID_(\d+)', jobs_list) if j == ID]
            if len(job) == 0:
                running = False
            else:
                print(f'Waiting for jobs to finish... {datetime.now():%Y-%m-%d %H:%M:%S}')
                time.sleep(sleep)
        print(f'Job ID_{ID} is finished')

    if load:
        os.makedirs(f'ssh_files/{ID}/', exist_ok=True)
        post_files = os.popen(f'ssh pegasus ls {job_folder}*{return_filetype}').read().split("\n")[:-1]

        if prejob_files is None:
            prejob_files = []
        import_files = [file for file in post_files if file not in prejob_files]

        for file in import_files:
            print('Copying {} to {}/ssh_files/{}{}/'.format(file, os.getcwd(), job_name, ID))
            os.system('scp pegasus:{} ssh_files/{}{}/{}'.format(file, job_name, ID, file.split('/')[-1]))
            image_display('ssh_files/{}{}/{}'.format(job_name, ID, file.split('/')[-1]))

    if check_IO_logs:
        logs = {'ErrorFile': '{}/*_logs_{}.stderr*'.format(job_folder, ID),
                'OutFile': '{}/*_logs_{}.stdout*'.format(job_folder, ID)
                }
        os.makedirs('logs/', exist_ok=True)
        for key, log in logs.items():
            os.system("scp 'pegasus:{}' 'logs/ID_{}_{}.txt'".format(log, ID, key))
            if os.path.isfile('logs/ID_{}_{}.txt'.format(ID, key)):
                print('logs/ID_{} {}:'.format(ID, key))
                with open('logs/ID_{}_{}.txt'.format(ID, key)) as file:
                    print(file.read())


def deeptools(regions, 
              signals, 
              matrix_name, 
              out_name, 
              pegasus_folder, 
              envelope='deeptools', 
              copy=False, 
              title='', 
              bps=(1500, 1500, 4000), 
              d_type='center', 
              scaled_names=('TSS', 'TES'), 
              make=('matrix', 'heatmap', 'heatmap_group', 'profile', 'profile_group'),
              missing_values_as_zero=True,
              heatmap_kmeans=0,
              save_sorted_regions='',
              sort_regions='descend'):
    '''
    Inputs
    ------
    regions: dictionary {'region_name':'/path/to/ssh/bedfile'}
    signals: dictionary {'signal_name':'/path/to/ssh/bigwigfile'}
    matrix_name: string of matrix name or matrix to be named (before .matrix.gz)
    out_name: name for output file
    tite: plot title (optional)
    envelope: conda envelope
    bps: tuple of region width on either side of center or scaled.  center ignores last number.  default is (1500,1500,4000)
    type: 'center' or 'scaled'
    scaled_names: optional names for scaled start and end (default ('TSS','TES'))
    make: tuple of deeptool commands.  options: matrix, heatmap, heatmap_group, profile, profile_group
    copy: bool.  Copy region and signal files to peagasus
    missing_values_as_zero: True
    heatmap_kmeans: Default 0.  kmeans clusters (int)
    save_sorted_regions= '' (default: don't output) else filename for kmeans sorted region file
    sort_regions= default descend.  'keep', 'no', ascend.
    Returns
    -------
    string of commands for ssh_job
    '''
    pegasus_folder = pegasus_folder if pegasus_folder.endswith('/') else f'{pegasus_folder}/'
    os.system(f"ssh pegasus 'mkdir {pegasus_folder}'")

    make_lower = [x.lower() for x in make]

    if d_type.lower() == 'center':
        deepMat = 'reference-point --referencePoint center'
        deepHeat = "--refPointLabel 'Peak Center'"
        deepProf = "--refPointLabel 'Peak Center'"
    else:
        deepMat = f'scale-regions --regionBodyLength {bps[2]}'
        deepHeat = f'--startLabel {scaled_names[0]} --endLabel {scaled_names[1]}'
        deepProf = f'--startLabel {scaled_names[0]} --endLabel {scaled_names[1]}'

    cmd_list = ['module rm python share-rpms65', f'source activate {envelope}']

    if copy:
        print('Copying region files to pegasus...')
        for region in regions.values():
            if os.popen(f'''ssh pegasus "if [ -f {pegasus_folder}{region.split('/')[-1]}]; then echo 'True' ; fi"''').read() != 'True\n':
                print(f'Copying {region} to pegasus at {pegasus_folder}.')
                os.system(f"scp {region} pegasus:{pegasus_folder}")
            else:
                print(f'{region} found in {pegasus_folder}.')

        print('Copying signal files to pegasus...')
        for signal in signals.values():
            if os.popen(f'''ssh pegasus "if [ -f {pegasus_folder}/{signal.split('/')[-1]} ]; then echo 'True' ; fi"''').read() != 'True\n':
                print(f'Copying {signal} to {pegasus_folder}.')
                os.system(f"scp {signal} pegasus:{pegasus_folder}")

        pegasus_region_path = ' '.join([f"{pegasus_folder}{region_path.split('/')[-1]}" for region_path in regions.values()])
        pegasus_signal_path = ' '.join([f"{pegasus_folder}{signal_path.split('/')[-1]}" for signal_path in signals.values()])
    else:
        pegasus_region_path = ' '.join([f'{region_path}' for region_path in regions.values()])
        pegasus_signal_path = ' '.join([f'{signal_path}' for signal_path in signals.values()])

    if 'matrix' in make_lower:
        signal_name = ' '.join([f'''"{signal_name.replace('_', ' ')}"''' for signal_name in signals.keys()])
        computeMatrix = f"computeMatrix {deepMat} -a {bps[0]} -b {bps[1]} -p 4 -R {pegasus_region_path} -S {pegasus_signal_path} --samplesLabel {signal_name} -o {matrix_name}.matrix.gz"
        if missing_values_as_zero:
            computeMatrix += ' --missingDataAsZero'
        cmd_list.append(computeMatrix)

    if 'heatmap' in make_lower or 'heatmap_group' in make_lower:
        region_name = ' '.join([f'''"{region_name.replace('_', ' ')}"''' for region_name in regions.keys()])
        plotHeatmap_base = f"plotHeatmap -m {matrix_name}.matrix.gz --dpi 300 {deepHeat} --plotTitle '{title.replace('_',' ')}' --whatToShow 'heatmap and colorbar' --colorMap Reds"
        if sort_regions != 'descend':
            plotHeatmap_base += f' --sortRegions {sort_regions}'
        if heatmap_kmeans > 0:
            plotHeatmap_base += f' --kmeans {heatmap_kmeans}'
        else:
            plotHeatmap_base += f' --regionsLabel {region_name}'
        if save_sorted_regions != '':
            plotHeatmap_base += f' --outFileSortedRegions {save_sorted_regions}.txt'
        if 'heatmap' in make_lower:
            cmd_list.append(f"{plotHeatmap_base} -out {out_name}_heatmap.png")
        if 'heatmap_group' in make_lower:
            cmd_list.append(f"{plotHeatmap_base} -out {out_name}_heatmap_perGroup.png --perGroup")

    if 'profile' in make_lower or 'profile_group' in make_lower:
        region_name = ' '.join([f'''"{region_name.replace('_', ' ')}"''' for region_name in regions.keys()])
        plotProfile_base = f"plotProfile -m {matrix_name}.matrix.gz --dpi 300 {deepProf} --plotTitle '{title.replace('_',' ')}'"
        if heatmap_kmeans > 0:
            plotProfile_base += f' --kmeans {heatmap_kmeans}'
        else:
            plotProfile_base += f' --regionsLabel {region_name}'
        if save_sorted_regions != '':
            plotProfile_base += f' --outFileSortedRegions {save_sorted_regions}_profile.txt'
        if 'profile' in make_lower:
            cmd_list.append(f"{plotProfile_base} -out {out_name}_profile.png")
        if 'profile_group' in make_lower:
            cmd_list.append(f"{plotProfile_base} -out {out_name}_profile_perGroup.png --perGroup")

    return cmd_list


def order_cluster(dict_set, count_df, gene_column_name, title):
    '''
    Inputs
    ------
    dict_set: a dictary with a cluster name and a set of genes in that cluster for plotting (should be non-overlapping).
    df: a pandas dataframe with the normalized counts for each gene and samples (or average of samples) in row columns.
              should also contain a column with the gene name.
    gene_column_name: the pandas column specifying the gene name (used in the dict_set)
    title: title for the plot and for saving the file
    Returns
    ------
    (Ordered Index List, Ordered Count DataFrame, Clustermap)
    '''
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.cluster import hierarchy
    import matplotlib.patches as mpatches

    out_list = []
    df = count_df.copy()
    df['group'] = 'NA'

    for name, genes in dict_set.items():
        if len(genes) == 0:
            print(f'There are not genes in {name}. Skipping Group')
            continue
        reduced_df = df[df[gene_column_name].isin(genes)]
        linkage = hierarchy.linkage(reduced_df.drop(columns=[gene_column_name, 'group']), method='ward', metric='euclidean')
        order = hierarchy.dendrogram(linkage, no_plot=True, color_threshold=-np.inf)['leaves']
        gene_list = reduced_df.iloc[order][gene_column_name].tolist()
        gene_index = df[df[gene_column_name].isin(gene_list)].index.tolist()
        out_list += gene_index

        gene_symbol = [gene.split('_')[-1] for gene in gene_list]
        with open(f'{name}_genes.txt', 'w') as file:
            for gene in gene_symbol:
                    file.write(f'{gene}\n')

        df.loc[gene_index, 'group'] = name

    ordered_df = df.loc[out_list]
    color_mapping = dict(zip([name for name, genes in dict_set.items() if len(genes) > 0], sns.hls_palette(len(df.group.unique()), s=.7)))
    row_colors = df.group.map(color_mapping)

    sns.set(context='notebook', font='Arial', palette='RdBu_r', style='white', rc={'figure.dpi': 300})
    clustermap = sns.clustermap(ordered_df.loc[out_list].drop(columns=[gene_column_name, 'group']),
                                z_score=0,
                                row_colors=row_colors,
                                row_cluster=False,
                                col_cluster=False,
                                cmap='RdBu_r',
                                yticklabels=False)
    clustermap.fig.suptitle(title)

    legend = [mpatches.Patch(color=color, label=label.replace('_', ' ')) for label, color in color_mapping.items() if label != 'NA']
    clustermap.ax_heatmap.legend(handles=legend, bbox_to_anchor=(-.1, .9, 0., .102))

    clustermap.savefig(f'{title.replace(" ","_")}.png', dpi=300)
    plt.close()
    image_display(f'{title.replace(" ","_")}.png')

    return out_list, ordered_df, clustermap


def ranked_ordered_cluster(dict_set, in_df, 
                           gene_column_name, 
                           dict_sort_col, 
                           title='ranked_ordered_cluster', 
                           group_name='Group',
                           figsize=None, 
                           ascending=False):
    '''
    Inputs
    ------
    dict_set: a dictary with a cluster name and a set of genes in that cluster for plotting.
    df: a pandas dataframe with the normalized counts for each gene and samples (or average of samples) in row columns.
              should also contain a column with the gene name.
    gene_column_name: the pandas column specifying the gene name (used in the dict_set)
    dict_sort_col: dictionary mapping cluster name with column to sort by in that cluster.
    group_name: name (string) of the clusters (ie. Group, or Lineage)
    title: title for the plot and for saving the file
    figsize: tuple of figsize or default none for autogeneration
    ascending: bool for sort order
    Returns
    ------
    (Ordered Count DataFrame, Clustermap)
    '''

    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats
    import matplotlib.patches as mpatches
    from dkfunctions import image_display

    out_dfs = []
    df = in_df.copy()
    df[group_name] = 'NA'
    df.index = df[gene_column_name]
    
    for name, genes in dict_set.items():
        reduced_df = df[df[gene_column_name].isin(genes)].copy()
        zscored = reduced_df.drop(columns=[gene_column_name, group_name]).T.apply(stats.zscore).T.copy()
        order = zscored.sort_values(by=dict_sort_col[name], ascending=ascending).index.tolist()
        gene_list = reduced_df.loc[order, gene_column_name].tolist()

        gene_symbol = [gene.split('_')[-1] for gene in gene_list]
        with open(f'{name}_genes.txt', 'w') as file:
            for gene in gene_symbol:
                    file.write(f'{gene}\n')

        reduced_df[group_name] = name
        reduced_df = reduced_df.loc[gene_list]
        out_dfs.append(reduced_df)

    ordered_df = pd.concat(out_dfs)
    
    groups = ordered_df[group_name].unique()
    color_mapping = dict(zip(groups, sns.color_palette("colorblind",len(groups))))
    row_colors = ordered_df[group_name].map(color_mapping).tolist()

    sns.set(context='paper', font='Arial', palette='pastel', style='white', rc={'figure.dpi': 300}, font_scale=.9)
    g = sns.clustermap(ordered_df.drop(columns=[gene_column_name, group_name]), 
                       z_score=0, 
                       row_colors=row_colors, 
                       row_cluster=False, 
                       col_cluster=False, 
                       cmap='RdBu_r', 
                       yticklabels=True,
                       figsize=figsize)
    
    g.fig.suptitle(title)
    
    legend = [mpatches.Patch(color=color, label=label.replace('_', ' ')) for label, color in color_mapping.items() if label != 'NA']
    g.ax_heatmap.legend(handles=legend, bbox_to_anchor=(-.1, .9, 0., .102),fontsize='large')
    g.savefig(f'{title.replace(" ","_")}.png', dpi=300)
    g.savefig(f'{title.replace(" ","_")}.svg')

    plt.close()
    image_display(f'{title.replace(" ","_")}.png')

    return ordered_df, g


def gsea_barplot(out_dir, pos_file, neg_file, gmt_name, max_number=20):
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
    import matplotlib.pyplot as plt
    import seaborn as sns

    out_dir = out_dir if out_dir.endswith('/') else '{}/'.format(out_dir)
    out_dir = '' if out_dir == '/' else out_dir
    os.makedirs(out_dir, exist_ok=True)
    pos = pd.read_table(pos_file).head(max_number) if os.path.isfile(pos_file) else pd.DataFrame(columns=['FDR q-val'])
    pos[gmt_name] = [' '.join(name.split('_')[1:]) for name in pos.NAME.tolist()]
    neg = pd.read_table(neg_file).head(max_number) if os.path.isfile(neg_file) else pd.DataFrame(columns=['FDR q-val'])
    neg[gmt_name] = [' '.join(name.split('_')[1:]) for name in neg.NAME.tolist()]

    sns.set(context='paper', font='Arial', font_scale=.9, style='white', rc={'figure.dpi': 300, 'figure.figsize': (8, 6)})
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2)
    fig.suptitle('{} GSEA enrichment\n(q<0.05, max {})'.format(gmt_name, max_number))

    if len(pos[pos['FDR q-val'] < 0.05]) > 0:
        UP = sns.barplot(data=pos[pos['FDR q-val'] < 0.05], x='NES', y=gmt_name, color='firebrick', ax=ax1)
        UP.set_title('Positive Enrichment')
        sns.despine()

    if len(neg[neg['FDR q-val'] < 0.05]) > 0:
        DN = sns.barplot(data=neg[neg['FDR q-val'] < 0.05], x='NES', y=gmt_name, color='steelblue', ax=ax2)
        DN.set_title('Negative Enrichment')
        sns.despine()

    try:
        plt.tight_layout(h_pad=1, w_pad=1)
    except ValueError:
        pass

    plt.subplots_adjust(top=0.88)
    file = f'{out_dir}{gmt_name}_GSEA_NES_plot.png'
    fig.savefig(file, dpi=300)
    plt.close()
    image_display(file)

    return file


def hinton(df, filename, folder, max_weight=None):
    """Draw Hinton diagram for visualizing a weight matrix."""

    import matplotlib.pyplot as plt
    import seaborn as sns

    folder = folder if folder.endswith('/') else f'{folder}/'
    folder = f'{os.getcwd()}/' if folder == '/' else folder

    sns.set(context='paper', rc={'figure.figsize': (8, 8), 'figure.dpi': 200})

    matrix = df.values
    plt.clf()
    plt.figure(figsize=(10, 10), dpi=200)
    ax = plt.gca()

    if not max_weight:
        max_weight = 2 ** np.ceil(np.log(np.abs(matrix).max()) / np.log(2))

    ax.patch.set_facecolor('white')
    ax.set_aspect('equal', 'box')
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.axis('off')

    for (x, y), w in np.ndenumerate(matrix):
        color = 'red' if w > 0 else 'blue'
        size = np.sqrt(np.abs(w) / max_weight)
        rect = plt.Rectangle([y - size / 2, x - size / 2], size, size,
                             facecolor=color, edgecolor=color)
        ax.add_patch(rect)

    fraction = len(df.index.tolist())
    increment = (.915 / fraction)
    y = 0.942
    for x in df.index.tolist():
        ax.annotate(x, xy=(-.15, y), xycoords='axes fraction')
        y -= increment
    ax.annotate("Components", xy=(.4, 0), xycoords='axes fraction', size=14)
    ax.autoscale_view()
    ax.annotate('Hinton Plot of Independent Components', xy=(.14, 1), xycoords='axes fraction', size=20)
    ax.invert_yaxis()
    ax.figure.savefig(f'{folder}{filename}.png')
    plt.close()
    image_display(f'{folder}{filename}.png')


def genomic_annotation_plots(dict_of_annotated_dfs, txdb_db,
                             filename='Genomic_Annotation_Plot',
                             title='',
                             bar_width=.75,
                             figsize=(10, 5),
                             order=['Promoter (<=1kb)',
                                    'Promoter (1-2kb)',
                                    'Promoter (2-3kb)',
                                    'Intron',
                                    'Exon',
                                    "3' UTR",
                                    "5' UTR",
                                    'Downstream (<1kb)',
                                    'Downstream (1-2kb)'
                                    'Downstream (2-3kb)',
                                    'Distal Intergenic'],
                             feature_col='annotation',
                             palette='colorblind',
                             plot_mode='fraction'
                             ):

    '''
    from chipseeker annotation output as df
    txdb_db = UCSC or Ensembl
    '''
    import matplotlib.pyplot as plt
    import seaborn as sns

    db = '(uc' if txdb_db == 'UCSC' else '(ENS'

    Anno_df = pd.DataFrame(index=order)

    for name, df in dict_of_annotated_dfs.items():
        df[feature_col] = [anno.replace(f' {db}', f'_{db}').split('_')[0] for anno in df[feature_col].tolist()]
        df_anno = df.groupby(feature_col).count().iloc[:, 0]
        if plot_mode.lower() == 'fraction':
            Anno_df[name] = df_anno / df_anno.sum()
        else:
            Anno_df[name] = df_anno

    Anno_df[Anno_df.isna()] = 0

    sns.set(style='white', font='Arial', font_scale=1.2)
    sns.set_palette(palette, n_colors=len(order))
    f = plt.figure(figsize=figsize)
    Anno_df.T.plot(kind='barh', stacked=True, ax=f.gca(), width=bar_width, lw=0.1)
    plt.title(title)
    plt.legend(loc=3, bbox_to_anchor=(1.0, 0))
    plt.xlabel('Fraction' if plot_mode.lower() == 'fraction' else 'Peak Number')
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300)
    plt.close()
    image_display(f'{filename}.png')


def extract_ENCODE_report_data(base_folder, report_type, out_folder='', histone=False, replicate=False):
    '''
    Inputs
    -----
    base_folder:  AQUAS results folder.  Will use subfolders for sample name and look for report in those subfolders.
    report_type: 'AQUAS' or 'cromwell'
    replicate: Whether the ChIPseq was performed as a repliate or not.
    Returns
    -----
    DataFrame of results
    '''

    tq = tq_type()

    if report_type.lower() not in ['aquas', 'cromwell']:
        raise ValueError('This function only extracts summary info from AQUAS or Cromwell generated qc reports.')

    base_folder = val_folder(base_folder)
    report_name = f'{base_folder}*/*report.html' if report_type.lower() == 'aquas' else f'{base_folder}*report.html'

    reports = glob.glob(report_name)
    out_folder = val_folder(out_folder)

    if replicate is True:
        raise AssertionError('Not set up for replicates yet.')

    results_df = pd.DataFrame(index=['Percent_mapped', 'Filtered_Uniquely_Mapped_Reads', 'Fraction_Duplicated', 'S_JS_Distance', 'PBC1', 'RSC', 'Overlap_Optimal_Peak_Number', 'FrIP_IDR', 'IDR_Peak_Number'])

    for file in tq(reports):
        name = re.findall(r'.*/(.*)_report.html', file)[0] if report_type.lower() == 'aquas' else re.findall(r'.*/(.*)_qc_report.html', file)[0]
        report = pd.read_html(file)
        series = pd.Series()
        series['Percent_mapped'] = report[1].iloc[7, 1] if report_type.lower() == 'aquas' else report[0].iloc[7, 1]
        series['Filtered_Uniquely_Mapped_Reads'] = report[2].iloc[5, 1] if report_type.lower() == 'aquas' else report[3].iloc[5, 1]
        series['Fraction_Duplicated'] = report[3].iloc[7, 1] if report_type.lower() == 'aquas' else report[1].iloc[7, 1]
        series['S_JS_Distance'] = report[4].iloc[7, 1] if report_type.lower() == 'aquas' else report[8].iloc[8, 1]
        series['PBC1'] = report[5].iloc[6, 1] if report_type.lower() == 'aquas' else report[2].iloc[6, 1]
        series['RSC'] = report[6].iloc[8, 1] if report_type.lower() == 'aquas' else report[5].iloc[9, 1]
        series['Overlap_Optimal_Peak_Number'] = report[10].iloc[4, 1] if report_type.lower() == 'aquas' else report[4].iloc[4, 1]

        if histone is False:
            series['FrIP_IDR'] = report[11].iloc[0, 1] if report_type.lower() == 'aquas' else report[7].iloc[1, 1]
            series['IDR_Peak_Number'] = report[12].iloc[4, 1] if report_type.lower() == 'aquas' else report[4].iloc[4, 2]
        results_df[name] = series

    for index in tq(results_df.index.tolist()):
        plot_col(results_df.loc[index], out=out_folder, title=f'{index}', ylabel=index.replace('_', ' '), plot_type=['violin', 'swarm'])

    return results_df


def meme_ssh(folder, fasta, bed, meme_db, out_name):
    folder = folder if folder.endswith('/') else f'{folder}/'
    out_fasta = f'{folder}{bed.split("/")[-1].replace(".bed",".fasta")}'

    meme_cmd = ['module rm python share-rpms65',
                'source activate motif',
                f'bedtools getfasta -fi {fasta} -bed {bed} -fo {out_fasta}',
                f'meme-chip -oc {out_name} -db {meme_db} -dna {out_fasta}'
                ]

    return ssh_job(meme_cmd, f'{out_name}_meme', folder, mem=3000)


def overlap_four(bed_dict, genome=None):
    '''
    Takes a dictionary of four pybedtools.BedTool objects.
    Merges all overlapping peaks for each bed into a master file.
    Intersects beds to merged master file.
    Performs annotations with ChIPseeker if genome is specified.
    Inputs
    ------
    bed_dict:  dictionary of BedTool files
    genome: 'hg38','hg19','mm10'
    Returns
    -------
    Returns a dictionary of dataframes from unique and overlap peaks.
    If genome is specified, includes a dictionary of annotated genes.
    '''
    from collections import OrderedDict
    import pickle

    names = list(bed_dict.keys())

    Folder = f'{os.getcwd()}/'
    subfolder = f"{names[0].replace(' ', '_')}-{ names[1].replace(' ', '_')}-{names[2].replace(' ', '_')}-{names[3].replace(' ', '_')}overlap/"

    out = f'{Folder}{subfolder}'
    os.makedirs(out, exist_ok=True)
    print(f'Output files are found in {out}')
    print(f'A: {names[0]}, B: {names[1]}, C: {names[2]}, D: {names[3]}')

    master = bed_dict[names[0]].cat(bed_dict[names[1]]).cat(bed_dict[names[2]]).sort().merge().cat(bed_dict[names[3]]).sort().merge()

    A = bed_dict[names[0]].sort().merge()
    B = bed_dict[names[1]].sort().merge()
    C = bed_dict[names[2]].sort().merge()
    D = bed_dict[names[3]].sort().merge()

    sorted_dict = OrderedDict({'master': master, 'A': A, 'B': B, 'C': C, 'D': D})
    sorted_dict['Abcd'] = (master + A - B - C - D)
    sorted_dict['aBcd'] = (master + B - A - C - D)
    sorted_dict['abCd'] = (master + C - A - B - D)
    sorted_dict['abcD'] = (master + D - A - B - C)

    sorted_dict['ABcd'] = (master + A + B - C - D)
    sorted_dict['AbCd'] = (master + A + C - B - D)
    sorted_dict['AbcD'] = (master + A + D - C - B)
    sorted_dict['aBCd'] = (master + B + C - A - D)
    sorted_dict['aBcD'] = (master + B + D - A - C)
    sorted_dict['abCD'] = (master + C + D - A - B)

    sorted_dict['ABCd'] = (master + A + B + C - D)
    sorted_dict['ABcD'] = (master + A + B + D - C)
    sorted_dict['AbCD'] = (master + A + C + D - B)
    sorted_dict['aBCD'] = (master + B + C + D - A)

    sorted_dict['ABCD'] = (master + A + B + C + D)

    labTup = tuple(key for key in sorted_dict.keys())
    lenTup = tuple(len(bed) for bed in sorted_dict.values())

    gener = (f'{lab}: {size}' for lab, size in zip(labTup, lenTup))
    for x in gener:
        print(x)

    for key, bed in sorted_dict.items():
        if len(bed) > 1:
            bed.to_dataframe().to_csv(f"{out}{key.replace(' ', '_')}-peaks-from-mergedPeaks.bed", header=None, index=None, sep="\t")

    if bool(genome):
        print('Annotating ovelapped peaks...')
        unikey = '{}'
        unianno = '{}_annotated'
        return_dict = annotate_peaks({unikey.format(key): bed.to_dataframe() for key, bed in sorted_dict.items() if len(bed) > 0}, out, genome=genome)

        gene_dict = {names[0]: return_dict[unianno.format('A')].SYMBOL.unique().tolist(),
                     names[1]: return_dict[unianno.format('B')].SYMBOL.unique().tolist(),
                     names[2]: return_dict[unianno.format('C')].SYMBOL.unique().tolist(),
                     names[3]: return_dict[unianno.format('D')].SYMBOL.unique().tolist()
                     }

        for name, gene_list in gene_dict.items():
            with open(f'{out}{name}_all_peaks_annotated.txt', 'w') as fp:
                for gene in gene_list:
                    fp.write(f'{gene}\n')

        with open(f'{out}Overlap_annotated_results.pkl', 'wb') as fp:
            pickle.dump(return_dict, fp)

    return sorted_dict if genome is None else {**sorted_dict, **return_dict}


def extract_clustermap_clusters(clustermap, num_of_clusters):

    '''
    Input a seaborn clustermap and number of clusters to id
    Returns an array of labelled clusters based on the original dataframe used to generate the clustermap.
    Usage: df['cluster'] = extract_clutsermap_clusters(clustermap, 2)
    '''
    from scipy.cluster.hierarchy import fcluster

    return fcluster(clustermap.dendrogram_row.linkage, num_of_clusters, criterion='maxclust')


def generate_ROSE_gffs(bed_file, name):
    '''
    takes a bed_file and returns a gff in the format needed for the ROSE algorithm.
    '''
    if type(bed_file) is BedTool:
        bed_file = bed_file.to_dataframe()
    elif type(bed_file) is not pd.DataFrame:
        IOError('Input must be either a pybedtools object or a pandas dataframe')

    gff = pd.DataFrame({'chr': bed_file.chrom,
                        'number': bed_file.index.tolist(),
                        'b1': '.',
                        'start': bed_file.start,
                        'end': bed_file.end,
                        'b2': '.',
                        'b3': '.',
                        'b4': '.',
                        'id': bed_file.index.tolist()
                        },
                       index=bed_file.index)

    gff.to_csv(f'{name}_enhancer.gff', header=False, index=False, sep="\t")
    return gff


def active_enhancer_determination(H3K4me1_bw, H3K4me3_bw, H3K4me1_bed, H3K27ac_bed, name, TSS_bed=None, gtf_file=None, TSS_region=(2500, 2500), chrom_sizes=None):
    '''
    Input
    ----------
    TSS_bed: location of a premade TSS window bed file.  If None, gtf_file, chrom.sizes, and TSS_region must be specified.
    gtf_file:  if no TSS_bed is provided, location of gtf_file must be included here.
    TSS_region: if no TSS_bed is provided, specify the distances before and after TSS.
    chrom: file of chromosome sizes (no header in file)
    H3K4me1_bw: location of H3K4me1 bigwig file (should be normalized ie. FC over input)
    H3K4me3_bw: location of H3K4me3 bigwig file (should be normalized ie. FC over input)
    H3K4me1_bed: peakfile for H3K4me1
    H3K27ac_bed: peakfile for H2K27ac
    name: name of the output sample
    Output
    -----------
    Saves enhancer bed file to disk as well as TSS window bed file if TSS_bed is not provided
    returns an enhancer bedfile
    '''

    import gtfparse
    import pyBigWig

    if TSS_bed is None:
        gtf = gtfparse.read_gtf(gtf_file)

        TSS = gtf[gtf.feature == 'gene'][['seqname', 'start', 'end', 'strand']]
        TSS['TSS'] = TSS[['start', 'end', 'strand']].apply(lambda x: x[0] if x[2] == '+' else x[1], axis=1)

        TSS_slop = pd.DataFrame({'chr': TSS.seqname,
                                 'start': TSS.TSS - TSS_region[0],
                                 'end': TSS.TSS + TSS_region[1]
                                 },
                                index=TSS.index)

        low_index = TSS_slop[TSS_slop.start < 0].index
        TSS_slop.loc[low_index, 'start'] = 0

        chrom = pd.read_csv(chrom_sizes, header=None, index_col=0, sep="\t")[1].to_dict()

        TSS_max = TSS_slop.groupby('chr').end.max().to_dict()

        problem_chroms = []
        for key, value in TSS_max.items():
            if value > chrom[key]:
                problem_chroms.append(key)

        for p_chrom in problem_chroms:
            idx = TSS_slop[TSS_slop.chr == p_chrom].index
            TSS_slop.loc[idx, 'end'] = TSS_slop.loc[idx, 'end'].apply(lambda x: x if x < chrom[p_chrom] else chrom[p_chrom])

        TSS_slop.to_csv(f'TSS_{TSS_region[0]}-{TSS_region[1]}.bed', index=False, header=False, sep="\t")
        TSS = BedTool.from_dataframe(TSS_slop)

    else:
        TSS = BedTool(TSS_bed)

    K1_bw = pyBigWig.open(H3K4me1_bw)
    K3_bw = pyBigWig.open(H3K4me3_bw)
    K1_bed = BedTool(H3K4me1_bed).to_dataframe()
    K27_bed = BedTool(H3K27ac_bed)

    df = K1_bed
    df['K4me1_mean'] = df.iloc[:, 0:3].apply(lambda x: K1_bw.stats(x[0], x[1], x[2], type='mean')[0], axis=1)
    df['K4me3_mean'] = df.iloc[:, 0:3].apply(lambda x: K3_bw.stats(x[0], x[1], x[2], type='mean')[0], axis=1)

    K4me1H_K4me3L = BedTool.from_dataframe(df[df.K4me1_mean > (df.K4me3_mean * 2)])

    enhancers = K27_bed + K4me1H_K4me3L - TSS

    print(f"Number of peaks with 2x more K4me1 than K4me3: {len(K4me1H_K4me3L)}")
    print(f"Number of enhancers using -{len(TSS_region[0])} to +{len(TSS_region[1])}  window: {len(enhancers)}")

    enhancers.saveas(f'{name}_Enhancers.bed')
    return enhancers


def boxplot_significance(x, y, data, type_test=None, __init__set=None):

    import matplotlib.pyplot as plt
    from decimal import Decimal
    from scipy.stats import ttest_ind, ks_2samp

    group_amount = data[x].unique()
    # min_value = np.min(data.groupby(x)[y].min().values)
    iqr = data.groupby(x)[y].describe()['75%'] - data.groupby(x)[y].describe()['25%']
    min_q1 = data.groupby(x)[y].describe()['25%'] - (1.5 * iqr)
    max_q3 = data.groupby(x)[y].describe()['75%'] + (1.5 * iqr)
    group_to_stats = data.groupby(x)

    n_min_q1 = np.min(min_q1.values)
    n_max_q3 = np.max(max_q3.values)

    perncentage_add = (n_max_q3 - n_min_q1) * 0.1
    print(perncentage_add, 'percent')
    print(group_amount)
    groups = {}
    if __init__set is not None:
        group_number = {g: (-0.25 + (0.5 * n)) + __init__set for n, g in enumerate(group_amount)}
    else:
        group_number = {g: n for n, g in enumerate(group_amount)}

    stack_lines = 1
    for g in group_amount:
        for i in group_amount:
            if g != i:
                if g + ':' + i not in groups and i + ':' + g not in groups:

                    print(g, i)
                    groups[g + ':' + i] = 0
                    n_g, n_i = group_number[g], group_number[i]

                    stack_value_plot = (perncentage_add) * stack_lines
                    print(stack_value_plot)
                    print('min_q1', n_min_q1)
                    plt.plot([n_g, n_i], [n_min_q1 - stack_value_plot, n_min_q1 - stack_value_plot], color='black')  # add the percentage

                    g_wo_na = group_to_stats.get_group(g)[y].dropna()
                    i_wo_na = group_to_stats.get_group(i)[y].dropna()

                    if type_test == 'ttest':
                        print(g_wo_na.median())
                        print(i_wo_na.median())
                        print('=' * 20)
                        test = ttest_ind(g_wo_na, i_wo_na)[1]
                    elif type_test == 'ks':

                        test = ks_2samp(g_wo_na, i_wo_na)[1]
                    else:
                        print('Please, enter  "ttest" or "ks" ')

                    print(n_i, n_g, test)

                    if test < 0.05:
                        sig_color = 'red'
                    else:
                        sig_color = 'black'

                    plt.text(((n_i - n_g) / 4) + n_g, n_min_q1 - (stack_value_plot * 0.9), '%.2E' % Decimal(test), color=sig_color)

                    stack_lines += 1

    plt.ylim(n_min_q1 - (perncentage_add * stack_lines), n_max_q3 + (perncentage_add))


def boxplot_significance_hue(x, y, hue, data, type_test=None):
    '''
    ################################################USAGE EXAMPLES##############################################
    # x,y boxplot
    sns.boxplot(x='rna_ratio_shift_groups',
                y='crt_mean',
                data=df_short_long_remove_inf, fliersize=0)
    boxplot_significance(x='rna_ratio_shift_groups',
                         y='crt_mean',
                         data=df_short_long_remove_inf, type_test='ttest')
    # x,y and  hue  boxplot
    sns.boxplot(x='rna_ratio_shift_groups',
                y='pause_index',
                hue='group',   # hue cagegory
                data=df_concat_pausebox,
                fliersize=0 )
    boxplot_significance_hue(x='rna_ratio_shift_groups', y='pause_index', hue='group', data=df_concat_pausebox, type_test='ttest' )
    '''

    data_hue = data.groupby(x)

    for e in enumerate(data_hue):
        k, d_h = e[1]
        n_loop = e[0]
        boxplot_significance(hue, y, d_h, type_test=type_test, __init__set=n_loop)


def logrank_km_hazard_ratio(km_df, time_col='months', censor_col='censor', group_col='group', reference=0):
    '''
    Determine the Hazard Ratio of a Kaplan Meier plot between to series in lifelines censor format

    Inputs
    ------
    km_df: dataframe of all samples with a column for the time, censor status, and which group (of two)
    time_col: column name for the time
    censor_col: column name for the censor status
    group_col: column name for the group
    reference: sample reference for HR (0 or 1)

    Outputs
    ------
    (Hazard Ratio, (Lower 95% CI, Upper 95% CI), Tau_df)

    '''
    from math import sqrt

    # Generate column of n_j (total number alive at each datapoint: 'cominbed_n_max')
    km_df.sort_values(by=time_col, inplace=True)
    km_df['combined_n'] = [x for x in range(len(km_df) - 1, -1, -1)]
    km_df = km_df.join(km_df.groupby(time_col).combined_n.max(), on=time_col, rsuffix='_max')

    # Generate column of n1_j and n2_j (total number of group1 and group2 alive at each datapoint)
    group1 = km_df[group_col].unique().tolist()[reference]
    group1_number = len(km_df[km_df[group_col] == group1])
    group1_n = []

    for x in km_df[group_col].tolist():
        in_group = 1 if x == group1 else 0
        group1_number -= in_group
        group1_n.append(group1_number)

    km_df['group1_n'] = group1_n
    km_df['group2_n'] = km_df.combined_n - km_df.group1_n

    # Make a new column with the number of each group alive at entering that timpoint (max of a groupby time)
    km_df = km_df.join(km_df.groupby(time_col).group1_n.max(), on=time_col, rsuffix='_max')
    km_df = km_df.join(km_df.groupby(time_col).group2_n.max(), on=time_col, rsuffix='_max')

    # Generate a column of d_j with total number of deaths at each timepoint
    km_df = km_df.join(km_df.groupby(time_col)[censor_col].sum(), on=time_col, rsuffix='_total_deaths')

    # Subset a dataframe with only non-censored data.
    km_data_tau = km_df[km_df[censor_col] == 1].copy()

    # Determine the expected values of group1 and group2 at each timepoint (Each relevant only within its own group)
    km_data_tau['e1'] = km_data_tau.group1_n_max * (km_data_tau[f'{censor_col}_total_deaths'] / km_data_tau.combined_n_max)
    km_data_tau['e2'] = km_data_tau.group2_n_max * (km_data_tau[f'{censor_col}_total_deaths'] / km_data_tau.combined_n_max)

    # Total Summed Expected and Observed values
    gr1_exp = km_data_tau.groupby(group_col).e1.sum()[0]
    gr2_exp = km_data_tau.groupby(group_col).e2.sum()[1]
    gr1_obs, gr2_obs = km_data_tau.groupby(group_col)[censor_col].sum()

    # Hazard Ratio
    HR = (gr1_obs / gr1_exp) / (gr2_obs / gr2_exp)

    # Standard Error of HR
    SE = sqrt((1 / gr1_exp) + (1 / gr2_exp))

    # Confidence intervals
    L = np.log(HR)
    L95 = np.exp(L - (1.96 * SE))
    U95 = np.exp(L + (1.96 * SE))

    return (HR, (L95, U95), km_data_tau)


def scale_factor(spike1_tags, spike2_tags, IP1_tags, IP2_tags):
    '''
    rep1 should be the lowest of number of tags in all samples in the comparison.
    returns scaled ratio of IP2.
    '''
    correction_factor = int(spike2_tags) / int(spike1_tags)
    scaled_IP2 = int(IP2_tags) / correction_factor
    IP2_ratio = scaled_IP2 / int(IP1_tags)

    return IP2_ratio


def spike_in_plot(spike_df, sample_list, description, out_dir):
    '''
    df: spike-in df with spike-in reads counts in column 'spike_reads' and genome counts in 'genome_reads'.
    sample_list:  list of samples to compare

    '''
    import seaborn as sns

    df = spike_df.loc[sample_list].copy()
    replicates = True if df.Replicate.unique.tolist() > 1 else False

    base_idx = df.genome_reads.astype(int).idxmin()
    base_spike = df.loc[base_idx, 'spike_reads']
    base_genome = df.loc[base_idx, 'genome_reads']

    df['Tag Ratio'] = df[['spike_reads', 'genome_reads']].apply(lambda x: scale_factor(base_spike, x[0], base_genome, x[1]))

    sns.set(context='notebook', style='white', font='Arial')
    if replicates:
        sns.boxplot(x='Condition', y='Tag Ratio', data=df, palette='pastel')
        sns.swarmplot(x='Condition', y='Ratio', hue='Replicate', data=df, size=8)
    else:
        sns.barplot(x='Condition', y='Ratio', data=df)

    sns.mpl.pyplot.ylabel('Normalized Genome Read Ratio')
    sns.mpl.pyplot.title(description.replace('_', ' '))
    sns.despine()
    sns.mpl.pyplot.savefig(f'{out_dir}{description.replaice(" ", "_")}.png')
    sns.mpl.pyplot.savefig(f'{out_dir}{description.replaice(" ", "_")}.svg')


def rename(file, name):
    if os.path.isfile(file):
        os.rename(file, name)


def globber(file):
    g = glob.glob(file)
    return '' if g == [] else g[0]


def ENCODE_clean(encode_folder, out_folder):
    '''
    folder is the ENCODE3 folder from chrome_chip that includes all crowmwell processed chipseq data.
    '''

    encode_folder = val_folder(encode_folder)
    out_folder = val_folder(out_folder)

    samples = os.listdir(encode_folder)

    os.makedirs(f'{out_folder}bams', exist_ok=True)
    os.makedirs(f'{out_folder}sample_peaks', exist_ok=True)
    os.makedirs(f'{out_folder}nodup_bam', exist_ok=True)
    os.makedirs(f'{out_folder}rep_overlap_peaks', exist_ok=True)
    os.makedirs(f'{out_folder}rep_idr_peaks', exist_ok=True)
    os.makedirs(f'{out_folder}reports', exist_ok=True)
    os.makedirs(f'{out_folder}bws', exist_ok=True)
    os.makedirs(f'{out_folder}jsons', exist_ok=True)

    for sample in samples:
        bw = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-macs2/shard-0/execution/*fc.signal.bigwig')
        sample_peak = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-macs2/shard-0/execution/*bfilt.narrowPeak.gz')
        nodup_bam = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-filter/shard-0/execution/*.nodup.bam')
        nodup_bai = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-filter/shard-0/execution/*.nodup.bam.bai')
        nodup_input_bam = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-filter_ctl/shard-0/execution/*nodup.bam')
        nodup_input_bai = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-filter_ctl/shard-0/execution/*nodup.bam.bai')
        bam = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-bwa/shard-0/execution/*.bam')
        bai = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-bwa/shard-0/execution/*.bam.bai')
        input_bam = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-bwa_ctl/shard-0/execution/*.bam')
        input_bai = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-bwa_ctl/shard-0/execution/*.bam.bai')
        optimal_peak = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz')
        conservative_peak = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-reproducibility_overlap/execution/conservative_peak.narrowPeak.gz')
        idr_optimal_peak = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz')
        idr_conservative_peak = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-reproducibility_idr/execution/conservative_peak.narrowPeak.gz')
        qc_report = globber(f'{encode_folder}{sample}/cromwell-executions/chip/*/call-qc_report/execution/qc.html')
        json = globber(f'{encode_folder}{sample}/*.json')

        rename(bw, f'{out_folder}bws/{sample}.fc.signal.bigwig')
        rename(sample_peak, f'{out_folder}sample_peaks/{sample}.500K.bfilt.narrowPeak.gz')
        rename(nodup_bam, f'{out_folder}nodup_bam/{sample}.nodup.bam')
        rename(nodup_bai, f'{out_folder}nodup_bam/{sample}.nodup.bam.bai')
        rename(nodup_input_bam, f'{out_folder}nodup_bam/{sample}_background.nodup.bam')
        rename(nodup_input_bai, f'{out_folder}nodup_bam/{sample}_background.nodup.bam.bai')
        rename(bam, f'{out_folder}bams/{sample}.bam')
        rename(bai, f'{out_folder}bams/{sample}.bam.bai')
        rename(input_bam, f'{out_folder}bams/{sample}_background.bam')
        rename(input_bai, f'{out_folder}bams/{sample}_background.bam.bai')
        rename(optimal_peak, f'{out_folder}rep_overlap_peaks/{sample}.optimal_peak.narrowPeak.gz')
        rename(conservative_peak, f'{out_folder}rep_overlap_peaks/{sample}.conservative_peak.narrowPeak.gz')
        rename(idr_optimal_peak, f'{out_folder}rep_idr_peaks/{sample}.idr.optimal_peak.narrowPeak.gz')
        rename(idr_conservative_peak, f'{out_folder}rep_idr_peaks/{sample}.idr.conservative_peak.narrowPeak.gz')
        rename(qc_report, f'{out_folder}reports/{sample}_qc_report.html')
        rename(json, f'{out_folder}jsons/{sample}_ENCODE3.json')


def normal_equation(X,y):
    '''
    OLS using the Normal Equation.
    inv(X.T.dot(X)).dot(X.T).dot(y)

    Inputs:
    X = array-like vector
    y = array-like vector

    Outputs:
    theta_paramaters(list), xspace(array), yspace(array)
    '''

    import numpy as np

    y = np.array(y)

    rows = int(np.size(X))

    #add the bias vector of ones
    X = np.hstack((np.ones((rows,1)),
                   np.reshape(x,(rows,1))
                  ))

    theta = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y)
    xspace = np.linspace(X[:,1].min(), X[:,1].max())
    yspace = np.array([theta[0] + theta[1]*x for x in xspace])

    return theta, xspace, yspace


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