#!/usr/local/python/bin/python
"""
---------------------------------------------
Create PRIZM data

Author: {author}
Contact: {email}
---------------------------------------------
"""

import subprocess
import glob
import string
import ConfigParser
import numpy as np
import pandas as pd
import re
import sys
import os
import operator
from operator import itemgetter
from collections import OrderedDict
from collections import defaultdict
import logging
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import multiprocessing
from multiprocessing import Process

__author__ = 'Hyukjung Kwon'
__email__ = "hjkwon@edgc.com"
__version__ = '0.9'

CFG_ReferenceMaleCountDir = ""
CFG_ReferenceFemaleCountDir = ""

logger = logging.getLogger(__name__)
# Specifies format of log
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter('%(asctime)s -%(lineno)s- [%(levelname) 9s] - %(message)s'))
logger.addHandler(ch)

# ==============================================================
def file_check(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file {0} does not exist!".format(arg))
    else:
        return str(arg)

# ==============================================================
def natural_keys(text):
    def atoi(text):
        return int(text) if text.isdigit() else text
    return [atoi(c) for c in re.split('(\d+)', text)]

# ==============================================================
def read_config(config_file):

    config = ConfigParser.RawConfigParser()
    config.read(config_file)

    # directory
    global CFG_ReferenceDir
    global CFG_ReferenceCountDir
    global CFG_ReferenceMaleCountDir
    global CFG_ReferenceFemaleCountDir
    global CFG_TrisomyControlDir
    global CFG_SampleRunDir

    CFG_ReferenceDir = config.get('directory', 'reference_dir')
    CFG_ReferenceCountDir = config.get('directory', 'reference_count_dir')
    CFG_ReferenceMaleCountDir = config.get('directory', 'reference_male_count_dir')
    CFG_ReferenceFemaleCountDir = config.get('directory', 'reference_female_count_dir')
    CFG_TrisomyControlDir = config.get('directory', 'trisomy_control_dir')
    CFG_SampleRunDir = config.get('directory', 'sample_run_count_dir')
    #CFG_ReferenceControlCountDir = config.get('directory', 'reference_control_count_dir')

    # data
    global CFG_CountExt
    global CFG_ZscoreExt
    #global CFG_Zscore10mbExt
    global CFG_CountSplitForHeatmap
    global CFG_CountSplitForGender
    global CFG_GenderList
    global CFG_MaleMean
    global CFG_MaleSd
    global CFG_FemaleMean
    global CFG_FemaleSd
    global CFG_Male10mbMean
    global CFG_Male10mbSd
    global CFG_Female10mbMean
    global CFG_Female10mbSd
    global CFG_TotalMean
    global CFG_TotalSd
    global CFG_Total10mbMean
    global CFG_Total10mbSd
    global CFG_TrisomyT13
    global CFG_TrisomyT18
    global CFG_TrisomyT21
    global CFG_TargetT13
    global CFG_TargetT18
    global CFG_TargetT21

    global CFG_NormalList
    global CFG_TrisomyList
    global CFG_DecisionAlgorithm

    CFG_CountExt = config.get('data', 'count_ext')
    CFG_ZscoreExt = config.get('data', 'zscore_ext')
    #CFG_Zscore10mbExt = config.get('data', 'zscore_10mb_ext')
    CFG_CountSplitForHeatmap = config.get('data', 'count_split_for_heatmap')
    CFG_CountSplitForGender = config.get('data', 'count_split_for_gender')
    CFG_GenderList = config.get('data', 'gender_list')
    CFG_MaleMean = config.get('data', 'male_mean')
    CFG_MaleSd = config.get('data', 'male_sd')
    CFG_FemaleMean = config.get('data', 'female_mean')
    CFG_FemaleSd = config.get('data', 'female_sd')
    CFG_Male10mbMean = config.get('data', 'male_10mb_mean')
    CFG_Male10mbSd = config.get('data', 'male_10mb_sd')
    CFG_Female10mbMean = config.get('data', 'female_10mb_mean')
    CFG_Female10mbSd = config.get('data', 'female_10mb_sd')
    CFG_TotalMean = config.get('data', 'total_mean')
    CFG_TotalSd = config.get('data', 'total_sd')
    CFG_Total10mbMean = config.get('data', 'total_10mb_mean')
    CFG_Total10mbSd = config.get('data', 'total_10mb_sd')
    CFG_NormalList = config.get('data', 'normal_list')
    CFG_TrisomyList = config.get('data', 'trisomy_list')
    CFG_TrisomyT13 = config.get('data', 'trisomy_control_t13')
    CFG_TrisomyT18 = config.get('data', 'trisomy_control_t18')
    CFG_TrisomyT21 = config.get('data', 'trisomy_control_t21')
    CFG_TargetT13 = config.get('data', 'target_t13')
    CFG_TargetT18 = config.get('data', 'target_t18')
    CFG_TargetT21 = config.get('data', 'target_t21')
    CFG_DecisionAlgorithm = config.get('src', 'multiz_core')

    print "\n======[C o n f i g u r a t i o n]======"
    print "reference_dir                : ", CFG_ReferenceDir
    print "reference_count_dir          : ", CFG_ReferenceCountDir
    print "reference_male_count_dir     : ", CFG_ReferenceMaleCountDir
    print "reference_female_count_dir   : ", CFG_ReferenceFemaleCountDir
    print "trisomy_control_dir          : ", CFG_TrisomyControlDir
    print ""
    print "count_ext                    : ", CFG_CountExt
    print "gender_list                  : ", CFG_GenderList
    print "trisomy_list                 : ", CFG_TrisomyList
    print "Decision algorithm           : ", CFG_DecisionAlgorithm
    print "=======================================\n"

# ==============================================================
def execute_cmd(command_list, shell_use):
    sp = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell_use)
    result, err = sp.communicate()
    if sp.returncode != 0:
        raise IOError(err)
        
# ==============================================================
def makeReference_10mb(file_list, mean_file, sd_file, mean_10mb_file, sd_10mb_file, autosomal):

    # 1. Reference by chromosome count
    dfs = {}
    for idx, filename in enumerate(file_list):
        count_data = pd.read_csv(filename, sep='\t')
        count_data = count_data[['space', 'start', 'reads']].rename(columns={'space':'chr', 'start':'bin', 'reads':'count'})
        count_data['bin'] = count_data['bin'].apply(lambda x: x/10000000)
        count_data_allsum = count_data.groupby(by=['chr'])['count'].sum()
        count_data_df = count_data_allsum.to_frame()
        sorted_index = sorted(count_data_df.index, key=natural_keys)
        count_data_df = count_data_df.loc[sorted_index]

        if autosomal == True:
            count_data_df = count_data_df.ix[0:22, :]

        count_sum = sum(count_data_df.ix[:,'count'])
        norm_data = [(count_data_df.ix[i,'count']/float(count_sum))*100 for i in range(len(count_data_df))]
        count_data_df.insert(1, 'ratio', norm_data)
        normalized_dict = {j : [count_data_df.ix[i,'ratio']/count_data_df.ix[j,'ratio'] for i in range(len(count_data_df))] for j in range(len(count_data_df))}
        dfs[idx] = pd.DataFrame.from_dict(normalized_dict)

    panel = pd.Panel(dfs)

    with open(mean_file, 'w') as f:
        mean_df = panel.mean(axis=0)
        logger.debug("Mean DF shape : %s " % str(mean_df.shape))
        mean_df.to_csv(f, sep='\t', index=False, header=False)
        logger.debug(mean_df.head())

    with open(sd_file, 'w') as f:
        sd_df = panel.std(axis=0)
        logger.debug("Sd DF shape : %s " % str(sd_df.shape)) 
        sd_df.to_csv(f, sep='\t', index=False, header=False)
        logger.debug(sd_df.head())

    # 2. Reference by 10mb bin count
    dfs = {}

    for idx, filename in enumerate(file_list):
        count_data = pd.read_csv(filename, sep='\t')
        count_data = count_data[['space', 'start', 'reads']].rename(columns={'space':'chr', 'start':'bin', 'reads':'count'})
        count_data['bin'] = count_data['bin'].apply(lambda x: x/10000000)
        count_data_allsum = count_data.groupby(by=['chr'])['count'].sum()
        count_data_df = count_data_allsum.to_frame()
        count_data_df.insert(0, 'bin', 0)
        sorted_index = sorted(count_data_df.index, key=natural_keys)

        count_data_df = count_data_df.drop(['chr9', 'chr13', 'chr16', 'chr18', 'chr21', 'chr22', 'chrX', 'chrY'])
        if autosomal == True:
            count_data_10mb = count_data.set_index('chr').loc[['chr9', 'chr13', 'chr16', 'chr18', 'chr21', 'chr22']]
        else:
            count_data_10mb = count_data.set_index('chr').loc[['chr9', 'chr13', 'chr16', 'chr18', 'chr21', 'chr22', 'chrX', 'chrY']]

        # merge two dataframes
        count_data_mix = pd.concat([count_data_df, count_data_10mb])
        if autosomal == True:
            count_data_mix = count_data_mix.loc[sorted_index[:-2]]
        else:
            count_data_mix = count_data_mix.loc[sorted_index]

        count_sum = sum(count_data_mix.ix[:,'count'])
        norm_data = [(count_data_mix.ix[i,'count']/float(count_sum))*100 for i in range(len(count_data_mix))]
        count_data_mix.insert(2, 'ratio', norm_data)

# this kind of dictionary comprehension is possible
        normalized_dict = {j : [count_data_mix.ix[i,'ratio']/count_data_mix.ix[j,'ratio'] if count_data_mix.ix[i, 'ratio'] != 0.0 and count_data_mix.ix[j, 'ratio'] != 0.0 else 0 for i in range(len(count_data_mix))] for j in range(len(count_data_mix))}

        dfs[idx] = pd.DataFrame.from_dict(normalized_dict)

    panel = pd.Panel(dfs)

    with open(mean_10mb_file, 'w') as f:
        mean_10mb_df = panel.mean(axis=0)
        logger.debug("Mean DF (10mb) shape : %s " % str(mean_10mb_df.shape)) 
        mean_10mb_df.to_csv(f, sep='\t', index=False, header=False)
        logger.debug(mean_10mb_df.head())

    with open(sd_10mb_file, 'w') as f:
        sd_10mb_df = panel.std(axis=0)
        logger.debug("Sd DF (10mb) shape : %s " % str(sd_10mb_df.shape)) 
        sd_10mb_df.to_csv(f, sep='\t', index=False, header=False)
        logger.debug(sd_10mb_df.head())

    return mean_df, sd_df, mean_10mb_df, sd_10mb_df

# ==============================================================
def calc_zscore_10mb(count_file, mean, sd, mean_10mb, sd_10mb):

    color_mapping_dict = {'chr9':'blue','chr13':'orange','chr16':'brown' ,'chr18':'purple','chr21':'red','chr22':'green','chrX': 'pink', 'chrY':'cyan'}

    logger.debug("calc_zscore_10mb run")
    # -------------------------------------
    # 1. Calculate zscore by chromosome
    # -------------------------------------
    count_data = pd.read_csv(count_file, sep='\t')
    count_data = count_data[['space', 'start', 'reads']].rename(columns={'space':'chr', 'start':'bin', 'reads':'count'})
    count_data['bin'] = count_data['bin'].apply(lambda x: x/10000000)
    count_data_allsum = count_data.groupby(by=['chr'])['count'].sum()
    count_data_df = count_data_allsum.to_frame()
    sorted_index = sorted(count_data_df.index, key=natural_keys)
    count_data_df = count_data_df.loc[sorted_index]

    count_sum = sum(count_data_df.ix[:,'count'])
    norm_data = [(count_data_df.ix[i,'count']/float(count_sum))*100 for i in range(len(count_data_df))]
    count_data_df.insert(1, 'ratio', norm_data)
    normalized_dict = {j : [count_data_df.ix[i,'ratio']/count_data_df.ix[j,'ratio'] for i in range(len(count_data_df))] for j in range(len(count_data_df))}
    normalized_df = pd.DataFrame(normalized_dict)

    z = {j : [(normalized_df.ix[i,j] - mean.ix[i,j])/sd.ix[i,j] if sd.ix[i,j] != 0.0 else 0.0 for i in range(len(count_data_df))] for j in range(len(count_data_df))}
    zscore = pd.DataFrame(z)
    #index_dict = {i:'chr%d'%(i+1) for i in range(22)}
    index_dict = {i:'chr%d'%(i+1) for i in range(22)}

    zscore = zscore.rename(index=index_dict, columns=index_dict)
    zscore.to_csv(count_file.split('10mb.')[0] + "zscore.txt", sep='\t')
    logger.debug(count_file.split('10mb.')[0] + "zscore.txt")
    index_dict.update({22:'chrX', 23:'chrY'})
    color_list = [color_mapping_dict.get(x, 'lightgray') for x in index_dict.values()]


    # -------------------------------------
    # 2. Calculate zscore by 10mb bin block
    # -------------------------------------

    # new column is at last column
    # drop 3 rows before merge
    count_data_df2 = count_data_allsum.to_frame()
    sorted_index = sorted(count_data_df2.index, key=natural_keys)
    count_data_df2 = count_data_df2.drop(['chr9', 'chr13', 'chr16', 'chr18', 'chr21', 'chr22', 'chrX', 'chrY'])
    count_data_df2.insert(0, 'bin', 0)
    count_data_10mb = count_data.set_index('chr').loc[['chr9', 'chr13', 'chr16', 'chr18', 'chr21', 'chr22', 'chrX', 'chrY']]

    # merge two dataframes
    count_data_mix = pd.concat([count_data_df2, count_data_10mb])
    count_data_mix = count_data_mix.loc[sorted_index]

    count_sum = sum(count_data_mix.ix[:,'count'])
    norm_data = [(count_data_mix.ix[i,'count']/float(count_sum))*100 for i in range(len(count_data_mix))]
    count_data_mix.insert(2, 'ratio', norm_data)

    # this kind of dictionary comprehension is possible
    normalizing = {j : [count_data_mix.ix[i,'ratio']/count_data_mix.ix[j,'ratio'] if count_data_mix.ix[i, 'ratio'] != 0.0 and count_data_mix.ix[j, 'ratio'] != 0.0 else 0 for i in range(len(count_data_mix))] for j in range(len(count_data_mix))}
    normalized_df = pd.DataFrame(normalizing)

    logger.debug(normalized_df.head())

    z = {j : [(normalized_df.ix[i,j] - mean_10mb.ix[i,j])/sd_10mb.ix[i,j] if sd_10mb.ix[i,j] != 0.0 else 0.0 for i in range(len(count_data_mix))] for j in range(len(count_data_mix))}
    zscore_10mb = pd.DataFrame(z)
    index_dict_10mb = {i:count_data_mix.index[i] for i in range(len(count_data_mix.index))}
    color_list_10mb = [color_mapping_dict.get(x, 'lightgray') for x in index_dict_10mb.values()]
    zscore_10mb = zscore_10mb.rename(index=index_dict_10mb, columns=index_dict_10mb)
    zscore_10mb.to_csv(count_file.split('10mb.')[0] + "zscore.10mb.txt", sep='\t')
    logger.debug(count_file.split('10mb.')[0] + "zscore.10mb.txt")
    logger.debug(zscore_10mb.head())

    row_count = len(count_data_mix)

    return zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count

# ==============================================================
def draw_heatmap(zscore, index_dict, color_list, ax, title):
    #ax = plt.gca()
    ax.set_xticks(np.arange(0, 24, 1))
    ax.set_yticks(np.arange(0, 24, 1))
    ax.set_title("[%s]" %(title), fontsize=18)

    ax.set_xticklabels(index_dict.values(), rotation=70, fontsize=12)
    ax.set_yticklabels(index_dict.values(), fontsize=12)

    # Minor ticks
    ax.set_xticks(np.arange(-.5, 24, 1), minor=True);
    ax.set_yticks(np.arange(-.5, 24, 1), minor=True);

    [t.set_color(i) for (i,t) in zip(color_list, ax.xaxis.get_ticklabels())]
    [t.set_color(i) for (i,t) in zip(color_list, ax.yaxis.get_ticklabels())]

    ax.imshow(zscore, interpolation='nearest', cmap=plt.get_cmap('RdBu_r'), vmin=-7, vmax=7)
    #ax.colorbar()
    ax.grid(which='minor', linestyle=':')
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')


# ==============================================================
def draw_heatmap_10mb(zscore, row_count, index_dict, color_list, ax):
    ax.set_xticks(np.arange(0, row_count, 1))
    ax.set_yticks(np.arange(0, row_count, 1))

    ax.set_xticklabels(index_dict.values(), rotation=70)
    ax.set_yticklabels(index_dict.values())

    # Minor ticks
    ax.set_xticks(np.arange(-.5, row_count, 1), minor=True);
    ax.set_yticks(np.arange(-.5, row_count, 1), minor=True);

    [t.set_color(i) for (i,t) in zip(color_list, ax.xaxis.get_ticklabels())]
    [t.set_color(i) for (i,t) in zip(color_list, ax.yaxis.get_ticklabels())]

    ax.imshow(zscore, interpolation='nearest', cmap=plt.get_cmap('RdBu_r'), vmin=-7, vmax=7)
    ax.grid(which='minor', linestyle=':', linewidth=0.3)

    ax.hlines(8-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(23-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(26-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(38-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(40-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(50-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(51-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(59-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(61-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(66-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(72-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.hlines(88-.5, -.5, row_count-.5, colors='black', linestyles="solid", linewidth= 1.0 )
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

# ==============================================================
def combine_dicts(a, b, op=operator.add):
    return dict(a.items() + b.items() + [(k, op(a[k], b[k])) for k in set(b) & set(a)])

# ==============================================================
def generateThresholdSet(Target_file, Trisomy_file, Trisomy_name):

    logger.info("make Thresholds for %s" %(Trisomy_name))

    if not os.path.exists(Target_file):
        # Male target file
        zscore_list = glob.glob(CFG_ReferenceMaleCountDir+"/*"+CFG_ZscoreExt)
        #logger.debug(zscore_list)
        extracted_data = []
        for zscore_file in zscore_list:
            filename = os.path.basename(zscore_file).split(CFG_CountSplitForGender)[0]
            #logger.debug(filename)
            with open(zscore_file, 'r') as f:
                for line in f:
                    a = line.split('\t')

                    if line.startswith('chr'+str(Trisomy_name[1:])):
                        extracted_data.append([filename, float(a[1]), float(a[2]), float(a[3]), float(a[4]), float(a[5]), float(a[6]), float(a[7]), float(a[8]), float(a[9]), float(a[10]), float(a[11]), float(a[12]), float(a[13]), float(a[14]), float(a[15]), float(a[16]), float(a[17]), float(a[18]), float(a[19]),float(a[20]), float(a[21]), float(a[22])])

        # Female target file
        zscore_list = glob.glob(CFG_ReferenceFemaleCountDir+"/*"+CFG_ZscoreExt)
        for zscore_file in zscore_list:
            filename = os.path.basename(zscore_file).split(CFG_CountSplitForGender)[0]
            #logger.debug(filename)
            with open(zscore_file, 'r') as f:
                for line in f:
                    a = line.split('\t')

                    if line.startswith('chr'+str(Trisomy_name[1:])):
                        extracted_data.append([filename, float(a[1]), float(a[2]), float(a[3]), float(a[4]), float(a[5]), float(a[6]), float(a[7]), float(a[8]), float(a[9]), float(a[10]), float(a[11]), float(a[12]), float(a[13]), float(a[14]), float(a[15]), float(a[16]), float(a[17]), float(a[18]), float(a[19]),float(a[20]), float(a[21]), float(a[22])])

        target_df = pd.DataFrame(extracted_data, columns=['sample name', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'])
        target_df.to_csv(Target_file, sep='\t') # for future use

        #logger.debug(target_df.head())

    else:
        target_df = pd.read_csv(Target_file, sep='\t', index_col=0)
        #logger.debug(target_df.head())
    
    if not os.path.exists(Trisomy_file):

        with open(CFG_TrisomyList) as f:
            trisomy_list = [tuple(x.strip().split()) for x in f.readlines()]
            trisomy_zipped_list = zip(*trisomy_list)
            #logger.debug(trisomy_zipped_list)

            file_list = [os.path.basename(x).split(CFG_CountSplitForGender)[0] for x in glob.glob(CFG_TrisomyControlDir+"/*"+CFG_ZscoreExt)]

            zscore_list = glob.glob(CFG_TrisomyControlDir+"/*"+CFG_ZscoreExt)
            #logger.debug(zscore_list)

            extracted_data = []
            for zscore_file in zscore_list:
                filename = os.path.basename(zscore_file).split(CFG_ZscoreExt)[0]
                filename_trisomysearch = os.path.basename(zscore_file).split(CFG_CountSplitForGender)[0]

                logger.debug(filename_trisomysearch)
                try:
                    if trisomy_zipped_list[1][trisomy_zipped_list[0].index(filename_trisomysearch)] == Trisomy_name:
                        logger.debug("%s --> %s" %(filename, Trisomy_name))

                        with open(zscore_file, 'r') as f:
                            for line in f:
                                a = line.split('\t')

                                if line.startswith('chr'+str(Trisomy_name[1:])):
                                    extracted_data.append([filename_trisomysearch, float(a[1]), float(a[2]), float(a[3]), float(a[4]), float(a[5]), float(a[6]), float(a[7]), float(a[8]), float(a[9]), float(a[10]), float(a[11]), float(a[12]), float(a[13]), float(a[14]), float(a[15]), float(a[16]), float(a[17]), float(a[18]), float(a[19]),float(a[20]), float(a[21]), float(a[22])])
                    else:
                        logger.debug("%s --> Not matched with '%s'" %(filename, Trisomy_name))

                except ValueError:
                    logger.debug("%s --> Not in Trisomy list" %(filename))
                    continue
            trisomy_df = pd.DataFrame(extracted_data, columns=['sample name', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'])
            trisomy_df.to_csv(Trisomy_file, sep='\t') # for future use
            #logger.debug(trisomy_df.head())

    else:
        trisomy_df = pd.read_csv(Trisomy_file, sep='\t', index_col=0)
        logger.debug(trisomy_df.head())

    # calculate & save thresholds
    min_trisomy = {chromosome : min(trisomy_df[chromosome]) for chromosome in trisomy_df.columns[1:]}
    max_normal = {chromosome : max(target_df[chromosome]) for chromosome in target_df.columns[1:]}

    diff_dic = combine_dicts(min_trisomy, max_normal, operator.sub)
    diff_list = sorted(diff_dic.items(), reverse=True, key=itemgetter(1))
    diff_list = [(chr, val) for (chr, val) in diff_list if val>-2.0]

    trisomy_list = [ min_trisomy[chr] for (chr, val) in diff_list ]
    normal_list = [ max_normal[chr] for (chr, val) in diff_list ]
    suspected_list = [np.average(e) for e in zip(trisomy_list, normal_list)]

    threshold_list = zip(trisomy_list, suspected_list, normal_list)
    threshold_dict = OrderedDict([(chr, threshold_list[i]) for i, (chr, val) in enumerate(diff_list)])

    plt.figure(dpi=150)
    plt.rcParams["figure.figsize"]= (14, 4)
    plt.rcParams['xtick.labelsize']=15
    plt.rcParams['ytick.labelsize']=15
    plt.rcParams['font.size']=8
    plt.rcParams['font.family']='sans-serif'

    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(22, 15, forward=True)

    plt.ylim(0,max(trisomy_list) + 1.0)
    plt.xlim(0,len(trisomy_list))

    x = np.arange(0, len(trisomy_list))
    plt.plot(x+0.5, trisomy_list, color='blue', linewidth=2.0)

    x = np.arange(0, len(suspected_list))
    plt.plot(x+0.5, suspected_list, color='orange', linewidth=2.0)

    x = np.arange(0, len(threshold_list))
    plt.plot(x+0.5, threshold_list, color='gray', linewidth=2.0)

    plt.xticks(x+0.5, [chr for (chr, val) in diff_list], rotation=60)
    plt.grid()
    plt.savefig("%s_thresholdset.png" %(Trisomy_name))
    plt.close()


    return threshold_dict

# ==============================================================
def run_zscore(count_list):

    name = multiprocessing.current_process().name
    logger.info("%s : Starting" % (name))

    with open(CFG_GenderList) as f:
        gender_list = [tuple(x.strip().split()) for x in f.readlines()]
        gender_zipped_list = zip(*gender_list)

    for count_file in count_list:

        filename_gendersearch = os.path.basename(count_file).split(CFG_CountSplitForGender)[0]
        # artificial sample_name
        # /data/syyun/NextSeq/Fastq/PRIZM/paired_10M2/control/artificial/2.98_ON17050420_S18.XYY.8p.3.10mb.txt
        #print count_file
        #filename_gendersearch = '_'.join(os.path.basename(count_file).split(CFG_CountSplitForGender)[1:]).split('.')[0]
        print filename_gendersearch
        try:
            if gender_zipped_list[1][gender_zipped_list[0].index(filename_gendersearch)] == 'XX':
                zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count= calc_zscore_10mb(count_file, female_mean_df, female_sd_df, female_10mb_mean_df, female_10mb_sd_df)
                logger.debug("%s --> XX" %(filename_gendersearch))
            else:
                zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count= calc_zscore_10mb(count_file, male_mean_df, male_sd_df, male_10mb_mean_df, male_10mb_sd_df)
                logger.debug("%s --> XY" %(filename_gendersearch))

            if args.draw_fig == 'Y':
                ax_heatmap = plt.subplot(gs[0:2,0:2])
                draw_heatmap(zscore, index_dict, color_list, ax_heatmap, filename_gendersearch)

                ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
                draw_heatmap_10mb(zscore_10mb, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)

                plt.tight_layout()
                filename_heatmap = os.path.basename(count_file).split(CFG_CountSplitForHeatmap)[0]
                plt.savefig("%s/%s.heatmap.png" %(CFG_SampleRunDir, filename_heatmap))
                logger.info("%s.heatmap.png created... OK" %(filename_heatmap))

        except ValueError:
            print "%s --> Not in the list" %(filename_gendersearch)
            continue

    logger.info("%s : Completed" % (name))


# ==============================================================
def run_algorithm(count_list):

    name = multiprocessing.current_process().name
    logger.info("%s : Starting" % (name))

    with open(CFG_GenderList) as f:
        gender_list = [tuple(x.strip().split()) for x in f.readlines()]
        gender_zipped_list = zip(*gender_list)

    for count_file in count_list:

        #filename_gendersearch = os.path.basename(count_file).split(CFG_CountSplitForGender)[0]
        #filename_gendersearch = '_'.join(os.path.basename(count_file).split(CFG_CountSplitForGender)[1:]).split('.')[0]
        filename_gendersearch = os.path.basename(count_file).split(CFG_CountSplitForGender)[0]
        #XXY = '_'.join(os.path.basename(count_file).split(CFG_CountSplitForGender)[1:]).split('.')[1]
        #logger.debug(XXY)
        try:
            # female
            #if gender_zipped_list[1][gender_zipped_list[0].index(filename_gendersearch)] == 'XX' and XXY != 'XXY':
            if gender_zipped_list[1][gender_zipped_list[0].index(filename_gendersearch)] == 'XX':
                logger.debug("%s --> XX" %(filename_gendersearch))
                command = "python " + CFG_DecisionAlgorithm + " -c10 " + count_file + " -m " + CFG_FemaleMean + " -s " + CFG_FemaleSd + " -m10 " + CFG_Female10mbMean + " -s10 " + CFG_Female10mbSd + " -single N -q 3.0"
                logger.debug(command)
                os.system(command)

            # male
            else:
                logger.debug("%s --> XY" %(filename_gendersearch))
                command = "python " + CFG_DecisionAlgorithm + " -c10 " + count_file + " -m " + CFG_MaleMean + " -s " + CFG_MaleSd + " -m10 " + CFG_Male10mbMean + " -s10 " + CFG_Male10mbSd + " -single N -q 3.0"
                logger.debug(command)
                os.system(command)

        except ValueError:
            print "%s --> Not in the list" %(filename_gendersearch)
            continue

    logger.info("%s : Completed" % (name))

# ==============================================================
def run_algorithm_noGender(count_list):

    name = multiprocessing.current_process().name
    logger.info("%s : Starting" % (name))
    for count_file in count_list:

        # female
        command = "python " + CFG_DecisionAlgorithm + " -c10 " + count_file + " -m " + CFG_FemaleMean + " -s " + CFG_FemaleSd + " -m10 " + CFG_Female10mbMean + " -s10 " + CFG_Female10mbSd + " -single N -q 3.0"
        logger.debug(command)
        os.system(command)

        '''
        # male
        else:
            logger.debug("%s --> XY" %(filename_gendersearch))
            command = "python " + CFG_DecisionAlgorithm + " -c10 " + count_file + " -m " + CFG_MaleMean + " -s " + CFG_MaleSd + " -m10 " + CFG_Male10mbMean + " -s10 " + CFG_Male10mbSd + " -lite Y -single N -q 3.0"
            logger.debug(command)
            os.system(command)
        '''

    logger.info("%s : Completed" % (name))

# ==============================================================
def run_heatmap(zscore_list):

    name = multiprocessing.current_process().name
    logger.info("%s : Starting" % (name))

    index_dict = {0: 'chr1', 1: 'chr2', 2: 'chr3', 3: 'chr4', 4: 'chr5', 5: 'chr6', 6: 'chr7', 7: 'chr8', 8: 'chr9', 9: 'chr10', 10: 'chr11', 11: 'chr12', 12: 'chr13', 13: 'chr14', 14: 'chr15', 15: 'chr16', 16: 'chr17', 17: 'chr18', 18: 'chr19', 19: 'chr20', 20: 'chr21', 21: 'chr22', 22: 'chrX', 23: 'chrY'}
    index_dict_10mb = {0: 'chr1', 1: 'chr2', 2: 'chr3', 3: 'chr4', 4: 'chr5', 5: 'chr6', 6: 'chr7', 7: 'chr8', 8: 'chr9', 9: 'chr9', 10: 'chr9', 11: 'chr9', 12: 'chr9', 13: 'chr9', 14: 'chr9', 15: 'chr9', 16: 'chr9', 17: 'chr9', 18: 'chr9', 19: 'chr9', 20: 'chr9', 21: 'chr9', 22: 'chr9', 23: 'chr10', 24: 'chr11', 25: 'chr12', 26: 'chr13', 27: 'chr13', 28: 'chr13', 29: 'chr13', 30: 'chr13', 31: 'chr13', 32: 'chr13', 33: 'chr13', 34: 'chr13',35: 'chr13', 36: 'chr13', 37: 'chr13', 38: 'chr14', 39: 'chr15', 40: 'chr16', 41: 'chr16', 42: 'chr16', 43: 'chr16', 44: 'chr16', 45: 'chr16', 46: 'chr16', 47: 'chr16', 48: 'chr16', 49: 'chr16', 50: 'chr17', 51: 'chr18', 52: 'chr18', 53: 'chr18', 54: 'chr18', 55: 'chr18', 56: 'chr18', 57: 'chr18', 58: 'chr18', 59: 'chr19', 60: 'chr20', 61: 'chr21', 62: 'chr21', 63: 'chr21', 64: 'chr21', 65: 'chr21', 66: 'chr22', 67: 'chr22', 68: 'chr22', 69: 'chr22', 70: 'chr22', 71: 'chr22', 72: 'chrX', 73: 'chrX', 74: 'chrX', 75: 'chrX', 76: 'chrX', 77: 'chrX', 78: 'chrX', 79: 'chrX', 80: 'chrX', 81: 'chrX', 82: 'chrX', 83: 'chrX', 84: 'chrX', 85: 'chrX', 86: 'chrX', 87: 'chrX', 88: 'chrY', 89: 'chrY', 90: 'chrY', 91: 'chrY', 92: 'chrY', 93: 'chrY'}

    color_mapping_dict = {'chr9':'blue','chr13':'orange','chr16':'brown' ,'chr18':'purple','chr21':'red','chr22':'green','chrX': 'pink', 'chrY':'cyan'}
    color_list = [color_mapping_dict.get(x, 'lightgray') for x in index_dict.values()]
    color_list_10mb = [color_mapping_dict.get(x, 'lightgray') for x in index_dict_10mb.values()]
    row_count = 94

    for zscore_file in zscore_list:
        ax_heatmap = plt.subplot(gs[0:2,0:2])
        #print zscore_file
        zscore_df = pd.read_csv(zscore_file, sep='\t')
        file_name = os.path.basename(zscore_file).split(CFG_CountSplitForHeatmap)[0]
        draw_heatmap(zscore_df, index_dict, color_list, ax_heatmap, file_name)

        ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
        zscore_10mb_file = '.'.join(zscore_file.split('.')[:-1]) + ".txt"
        zscore_10mb_df = pd.read_csv(zscore_10mb_file, sep='\t')
        draw_heatmap_10mb(zscore_10mb_df, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)

        plt.tight_layout()
        plt.savefig("%s/%s.heatmap.png" %(CFG_SampleRunDir, file_name))
        logger.info("%s.heatmap.png created... OK" %(file_name))

    logger.info("%s : Completed" % (name))

# ==============================================================
def split_seq_by_num_of_size(iterable, size):
    it = iter(iterable)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size))

def split_seq_by_num_of_sublist(iterable, num_of_sublist):
    k, m = divmod(len(iterable), num_of_sublist)
    return (iterable[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(num_of_sublist))

# ==============================================================
def doParallel(func, count_list):
    proc = []

    for i in range(len(count_list)):

        p = Process(name=str(i+1), target=func, args=(count_list[i],))
        p.start()
        proc.append(p)

    for p in proc:
        p.join()

# ==============================================================
#   MAIN FUNCTION   
# ==============================================================
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__.format(author=__author__, email=__email__), formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False)

    pgroup = parser.add_argument_group("Input")
    pgroup.add_argument('config', metavar='config', type=lambda x:file_check(parser,x), help='config file (setup.ini)')

    ogroup = parser.add_argument_group("Options")
    ogroup.add_argument('-classify', dest='classify', metavar='Y|N',default='N', help='Classify by gender')
    ogroup.add_argument('-reference', dest='reference', metavar='Y|N',default='N', help='Make reference file')
    ogroup.add_argument('-zscore_normal', dest='zscore_n', metavar='Y|N',default='N', help='Make reference zscore')
    ogroup.add_argument('-zscore_trisomy', dest='zscore_t', metavar='Y|N',default='N', help='Make trisomy zscore')
    ogroup.add_argument('-threshold', dest='threshold', metavar='Y|N',default='N', help='Make threshold')
    ogroup.add_argument('-draw_fig', dest='draw_fig', metavar='Y|N',default='N', help='draw & save heatmap')
    ogroup.add_argument('-target_dir', dest='target_dir', metavar='Y|N',default='N', help='calc zscore & save heatmap')
    ogroup.add_argument('-target_file', dest='target_file', type=lambda x:file_check(parser,x), help='calc zscore & save heatmap for individual count file')
    ogroup.add_argument('-decision', dest='decision', metavar='Y|N',default='N', help='run decision algorithm')
    ogroup.add_argument('-decision_ng', dest='decision_ng', metavar='Y|N',default='N', help='run decision algorithm no gender info')
    ogroup.add_argument('-heatmap', dest='heatmap', metavar='Y|N',default='N', help='draw heatmap for zscore')
    ogroup.add_argument('-core', dest='core', default=1, type=int, help='number of core to use for target zscore')
    ogroup.add_argument('-v','--version', action='version', version='%(prog)s '+ __version__)
    ogroup.add_argument('-h','--help',action='help', help='show this help message and exit')
    ogroup.add_argument('--debug', dest='debug', action='store_true', default=False, help=argparse.SUPPRESS)

    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)    

    pd.set_option('display.max_rows', 1000)
    pd.set_option('display.max_columns', 500)


    read_config(args.config)

    # -------------------------------------------------------
    # 1. classify count files according to gender information
    # -------------------------------------------------------

    if args.classify == 'Y':
        if os.path.isdir(CFG_ReferenceMaleCountDir) == False:
            execute_cmd(["mkdir", CFG_ReferenceMaleCountDir], False)

        if os.path.isdir(CFG_ReferenceFemaleCountDir) == False:
            execute_cmd(["mkdir", CFG_ReferenceFemaleCountDir], False)

        with open(CFG_GenderList) as f:
            gender_list = [tuple(x.strip().split()) for x in f.readlines()]
            gender_zipped_list = zip(*gender_list)

            #count_list = glob.glob(CFG_ReferenceCountDir +  "/*" +CFG_CountExt)
            count_list = [os.path.basename(x).split(CFG_CountSplitForGender)[0] for x in glob.glob(CFG_ReferenceCountDir+"/*"+CFG_CountExt)]

            for idx, filename in enumerate(count_list):

                #real_filename = os.path.basename(filename).split(CFG_CountSplit)[0]
                try:
                    if gender_zipped_list[1][gender_zipped_list[0].index(filename)] == 'XX':
                        logger.debug("%s --> XX" %(filename))
                        command_string = "cp %s/%s* %s" % (CFG_ReferenceCountDir, filename, CFG_ReferenceFemaleCountDir)
                    else:
                        logger.debug("%s --> XY" %(filename))
                        command_string = "cp %s/%s* %s" % (CFG_ReferenceCountDir, filename, CFG_ReferenceMaleCountDir)

                    execute_cmd(command_string, True)

                except ValueError:
                    print "%s --> Not in the list" %(filename)
                    continue
            
            #if filename in (item[0] for item in gender_list):

    # -------------------------------------------------------
    # 2. make reference
    # -------------------------------------------------------
    if args.reference == 'Y':

        # 2-1. make overall reference
        logger.info("Generate overall reference with %s" %(CFG_ReferenceCountDir))
        file_list = glob.glob(CFG_ReferenceCountDir+"/*"+CFG_CountExt)
        total_mean_df, total_sd_df, total_10mb_mean_df, total_10mb_sd_df = makeReference_10mb(file_list, CFG_TotalMean, CFG_TotalSd, CFG_Total10mbMean, CFG_Total10mbSd, autosomal=True)

        logger.info("Generate female reference with %s" %(CFG_ReferenceFemaleCountDir))
        file_list = glob.glob(CFG_ReferenceFemaleCountDir+"/*"+CFG_CountExt)
        female_mean_df, female_sd_df, female_10mb_mean_df, female_10mb_sd_df = makeReference_10mb(file_list, CFG_FemaleMean, CFG_FemaleSd, CFG_Female10mbMean, CFG_Female10mbSd, autosomal=False)

        # 2-3 combine two dataframe
        logger.debug(female_mean_df.head())
        female_mean_df.update(total_mean_df)
        female_sd_df.update(total_sd_df)

        logger.debug(female_10mb_mean_df.head())
        female_10mb_mean_df.update(total_10mb_mean_df)
        female_10mb_sd_df.update(total_10mb_sd_df)

        # 2-4 overwrite csv files
        with open(CFG_FemaleMean, 'w') as f:
            female_mean_df.to_csv(f, sep='\t', index=False, header=False)
        with open(CFG_FemaleSd, 'w') as f:
            female_sd_df.to_csv(f, sep='\t', index=False, header=False)
        with open(CFG_Female10mbMean, 'w') as f:
            female_10mb_mean_df.to_csv(f, sep='\t', index=False, header=False)
        with open(CFG_Female10mbSd, 'w') as f:
            female_10mb_sd_df.to_csv(f, sep='\t', index=False, header=False)


        logger.debug(female_sd_df.head())
        logger.debug("Female Sd DF : %s" % str(female_sd_df.shape))
        logger.debug(female_10mb_mean_df.head())
        logger.debug("Female 10mb Sd DF : %s" % str(female_10mb_sd_df.shape))

    else:
        logger.info("Read reference files")
        male_mean_df = pd.read_csv(CFG_MaleMean, sep='\t', header=None)
        male_sd_df = pd.read_csv(CFG_MaleSd, sep='\t', header=None)
        male_10mb_mean_df = pd.read_csv(CFG_Male10mbMean, sep='\t', header=None)
        male_10mb_sd_df = pd.read_csv(CFG_Male10mbSd, sep='\t', header=None)

        female_mean_df = pd.read_csv(CFG_FemaleMean, sep='\t', header=None)
        female_sd_df = pd.read_csv(CFG_FemaleSd, sep='\t', header=None)
        female_10mb_mean_df = pd.read_csv(CFG_Female10mbMean, sep='\t', header=None)
        female_10mb_sd_df = pd.read_csv(CFG_Female10mbSd, sep='\t', header=None)

    # -------------------------------------------------------
    # 4. create normal & trisomy control z-score
    # -------------------------------------------------------
    plt.figure(dpi=150)
    plt.rcParams["figure.figsize"]= (14, 4)
    plt.rcParams['xtick.labelsize']=5
    plt.rcParams['ytick.labelsize']=5
    plt.rcParams['font.size']=8
    plt.rcParams['font.family']='sans-serif'


    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(22, 15, forward=True)

    gs = GridSpec(4, 2)

    # 4-1. Normal zscore
    if args.zscore_n == 'Y':
        logger.info("Generate male normal z-score with %s" %(CFG_ReferenceMaleCountDir))
        file_list = glob.glob(CFG_ReferenceMaleCountDir+"/*"+CFG_CountExt)

        logger.info(file_list)

        for count_file in file_list:
            zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count= calc_zscore_10mb(count_file, male_mean_df, male_sd_df, male_10mb_mean_df, male_10mb_sd_df)

            if args.draw_fig == 'Y':
                ax_heatmap = plt.subplot(gs[0:2,0:2])
                file_name = os.path.basename(count_file).split(CFG_CountSplitForHeatmap)[0]
                draw_heatmap(zscore, index_dict, color_list, ax_heatmap, file_name)

                logger.debug(zscore.head())

                ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
                draw_heatmap_10mb(zscore_10mb, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)
                logger.debug(zscore_10mb.head())

                plt.tight_layout()
                plt.savefig("%s/%s.heatmap.png" %(CFG_ReferenceMaleCountDir, file_name))
                logger.info("%s.heatmap.png created... OK" %(file_name))


        logger.info("Generate female normal z-score with %s" %(CFG_ReferenceFemaleCountDir))
        file_list = glob.glob(CFG_ReferenceFemaleCountDir+"/*"+CFG_CountExt)

        for count_file in file_list:
            zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count= calc_zscore_10mb(count_file, female_mean_df, female_sd_df, female_10mb_mean_df, female_10mb_sd_df)

            if args.draw_fig == 'Y':
                ax_heatmap = plt.subplot(gs[0:2,0:2])
                file_name = os.path.basename(count_file).split(CFG_CountSplitForHeatmap)[0]
                draw_heatmap(zscore, index_dict, color_list, ax_heatmap, file_name)

                ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
                draw_heatmap_10mb(zscore_10mb, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)

                plt.tight_layout()
                plt.savefig("%s/%s.heatmap.png" %(CFG_ReferenceFemaleCountDir, file_name))
                logger.info("%s.heatmap.png created... OK" %(file_name))

    # 4-2. Trisomy zscore
    if args.zscore_t == 'Y':
        logger.info("Generate trisomy z-score with %s" %(CFG_TrisomyControlDir))
        file_list = glob.glob(CFG_TrisomyControlDir+"/*"+CFG_CountExt)

        with open(CFG_GenderList) as f:
            gender_list = [tuple(x.strip().split()) for x in f.readlines()]
            gender_zipped_list = zip(*gender_list)

        for count_file in file_list:
            filename_gendersearch = os.path.basename(count_file).split(CFG_CountSplitForGender)[0]
            try:
                if gender_zipped_list[1][gender_zipped_list[0].index(filename_gendersearch)] == 'XX':
                    logger.info("%s : XX" %(filename_gendersearch))
                    zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count= calc_zscore_10mb(count_file, female_mean_df, female_sd_df, female_10mb_mean_df, female_10mb_sd_df)
                else:
                    logger.info("%s : XY" %(filename_gendersearch))
                    zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count= calc_zscore_10mb(count_file, male_mean_df, male_sd_df, male_10mb_mean_df, male_10mb_sd_df)

                if args.draw_fig == 'Y':
                    ax_heatmap = plt.subplot(gs[0:2,0:2])
                    file_name = os.path.basename(count_file).split(CFG_CountSplitForHeatmap)[0]
                    draw_heatmap(zscore, index_dict, color_list, ax_heatmap, file_name)

                    logger.debug(zscore.head())

                    ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
                    draw_heatmap_10mb(zscore_10mb, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)
                    logger.debug(zscore_10mb.head())

                    plt.tight_layout()
                    plt.savefig("%s/%s.heatmap.png" %(CFG_TrisomyControlDir, file_name))
                    logger.info("%s.heatmap.png created... OK" %(file_name))

            except ValueError:
                print "%s --> Not in the list" %(filename_gendersearch)
                continue

    # -------------------------------------------------------
    # 5. create thresholds
    # -------------------------------------------------------
    if args.threshold == 'Y':
        logger.info("make Thresholds for each Trisomy")
        logger.info("Normal reference (male) : %s" %(CFG_ReferenceMaleCountDir))
        logger.info("Normal reference (female) : %s" %(CFG_ReferenceFemaleCountDir))
        logger.info("Trisomy reference (male+female) : %s" %(CFG_TrisomyControlDir))

        T13_thresholdDic = generateThresholdSet(CFG_TargetT13, CFG_TrisomyT13, "T13")
        print T13_thresholdDic
        T18_thresholdDic = generateThresholdSet(CFG_TargetT18, CFG_TrisomyT18, "T18")
        print T18_thresholdDic
        T21_thresholdDic = generateThresholdSet(CFG_TargetT21, CFG_TrisomyT21, "T21")
        print T21_thresholdDic


    # -------------------------------------------------------
    # 6. calculate z-score & draw heatmap (directory)
    # -------------------------------------------------------
    if args.target_dir == 'Y':
        logger.info("Generate z-score for target directory %s" %(CFG_SampleRunDir))

        #with open(CFG_GenderList) as f:
        #    gender_list = [tuple(x.strip().split()) for x in f.readlines()]
        #    gender_zipped_list = zip(*gender_list)

        #count_list = glob.glob(CFG_ReferenceCountDir +  "/*" +CFG_CountExt)
        #file_list = [os.path.basename(x).split(CFG_CountSplitForHeatmap)[0] for x in glob.glob(CFG_SampleRunDir+"/*"+CFG_CountExt)]

        count_list = glob.glob(CFG_SampleRunDir+"/*"+CFG_CountExt)
        split_list = split_seq_by_num_of_sublist(count_list, args.core)
        doParallel(run_zscore, list(split_list))

    # -------------------------------------------------------
    # 6-1. calculate z-score & draw heatmap (target_file)
    # -------------------------------------------------------
    if args.target_file is not None:
        logger.info("Generate z-score for target file %s" %(args.target_file))

        with open(CFG_GenderList) as f:
            gender_list = [tuple(x.strip().split()) for x in f.readlines()]
            gender_zipped_list = zip(*gender_list)

        filename_gendersearch = os.path.basename(args.target_file).split(CFG_CountSplitForGender)[0]
        try:
            if gender_zipped_list[1][gender_zipped_list[0].index(filename_gendersearch)] == 'XX':
                logger.info("%s : XX" %(filename_gendersearch))
                zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count= calc_zscore_10mb(args.target_file, female_mean_df, female_sd_df, female_10mb_mean_df, female_10mb_sd_df)
            else:
                logger.info("%s : XY" %(filename_gendersearch))
                zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count= calc_zscore_10mb(args.target_file, male_mean_df, male_sd_df, male_10mb_mean_df, male_10mb_sd_df)

            if args.draw_fig == 'Y':
                ax_heatmap = plt.subplot(gs[0:2,0:2])
                file_name = os.path.basename(args.target_file).split(CFG_CountSplitForHeatmap)[0]
                draw_heatmap(zscore, index_dict, color_list, ax_heatmap, file_name)

                logger.debug(zscore.head())

                ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
                draw_heatmap_10mb(zscore_10mb, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)
                logger.debug(zscore_10mb.head())

                plt.tight_layout()
                plt.savefig("%s/%s.heatmap.png" %(CFG_ReferenceMaleCountDir, file_name))
                logger.info("%s.heatmap.png created... OK" %(file_name))

        except ValueError:
            print "%s --> Not in the list" %(filename_gendersearch)

    # -------------------------------------------------------
    # 7. run multiz core algorithm (for each directory)
    # -------------------------------------------------------

    if args.decision == 'Y':
        logger.info("run multiZ core algorithm for %s" %(CFG_SampleRunDir))

        #with open(CFG_GenderList) as f:
        #    gender_list = [tuple(x.strip().split()) for x in f.readlines()]
        #    gender_zipped_list = zip(*gender_list)

        #    file_list = [os.path.basename(x).split(CFG_CountSplitForHeatmap)[0] for x in glob.glob(CFG_SampleRunDir+"/*"+CFG_CountExt)]
        count_list = glob.glob(CFG_SampleRunDir+"/*"+CFG_CountExt)
        split_list = split_seq_by_num_of_sublist(count_list, args.core)
        doParallel(run_algorithm, list(split_list))

    # -------------------------------------------------------
    # 8. run multiz core algorithm (without Gender)
    # -------------------------------------------------------

    if args.decision_ng == 'Y':
        logger.info("run multiZ core algorithm for %s" %(CFG_SampleRunDir))

        file_list = [os.path.basename(x).split(CFG_CountSplitForHeatmap)[0] for x in glob.glob(CFG_SampleRunDir+"/*"+CFG_CountExt)]
        count_list = glob.glob(CFG_SampleRunDir+"/*"+CFG_CountExt)

        split_list = split_seq_by_num_of_sublist(count_list, args.core)
        doParallel(run_algorithm_noGender, list(split_list))

    # -------------------------------------------------------
    # 9. draw heatmap only
    # -------------------------------------------------------
    if args.heatmap == 'Y':
        logger.info("Draw heatmap for targetted directory %s" %(CFG_SampleRunDir))

        zscore_list = glob.glob(CFG_SampleRunDir+"/*"+CFG_ZscoreExt)
        split_list = split_seq_by_num_of_sublist(zscore_list, args.core)
        doParallel(run_heatmap, list(split_list))

