#!/usr/local/python/bin/python
"""
---------------------------------------------
Create PRIZM data (10mb all)
chromosome(autosome) X 10mb all chromosome

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
def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

# ==============================================================
def getNormalized10mbData(filename, autosomal, key_start=0):

    df = pd.read_csv(filename, sep='\t')
    df = df[['space', 'start', 'reads']].rename(columns={'space':'chr', 'start':'bin', 'reads':'count'})
    df['bin'] = df['bin'].apply(lambda x: (x-1)/10000000)
    df.bin = df.bin.astype(int)
    # HMMcopy.R makes some different order, so I need to sort them again
    df = df.set_index('chr').loc[['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']]
    df.reset_index(level=0, inplace=True)
    # sum function doesn't follow the original order, so I need to sort it too
    count_data_allsum = df.groupby(by=['chr'])['count'].sum()
    sum_df = count_data_allsum.to_frame()
    sorted_index = sorted(sum_df.index, key=natural_keys)
    sum_df = sum_df.loc[sorted_index]
    #
    # df : 10mb dataframe
    # sum_df : chromosom dataframe
    #
    sum_df = sum_df.iloc[0:22, :] 
    count_sum = sum_df['count'].sum()
    #norm_data = [(sum_df.ix[i,'count']/float(count_sum))*100 for i in range(len(sum_df))]
    # new pandas uses loc/iloc
    norm_data = [(sum_df.loc[sum_df.index[i],'count']/float(count_sum))*100 for i in range(len(sum_df))]
    sum_df.insert(1, 'ratio', norm_data)

    if autosomal == True:
        df = df.loc[~df['chr'].isin(['chrX','chrY'])]
    else:
        df = df.loc[df['chr'].isin(['chrX','chrY'])]

    norm_10mb_data = [(df.loc[df.index[i],'count']/float(count_sum))*100 for i in range(len(df))]
    df.insert(2, 'ratio', norm_10mb_data)

    #normalized_dict = {j+key_start : [sum_df.loc[sum_df.index[i],'ratio']/df.loc[df.index[j],'ratio'] if sum_df.loc[sum_df.index[i], 'ratio'] != 0.0 and df.loc[df.index[j], 'ratio'] != 0.0 else 0 for i in range(len(sum_df))] for j in range(len(df))}
    normalized_dict = {j+key_start : [df.loc[df.index[j],'ratio']/sum_df.loc[sum_df.index[i],'ratio'] if sum_df.loc[sum_df.index[i], 'ratio'] != 0.0 and df.loc[df.index[j], 'ratio'] != 0.0 else 0 for i in range(len(sum_df))] for j in range(len(df))}

    return normalized_dict


# ==============================================================
def makeReferenceData(file_list, autosomal):

    dfs = {}
    for idx, filename in enumerate(file_list):
        normalized_dict = getNormalized10mbData(filename, autosomal)
        dfs[idx] = pd.DataFrame.from_dict(normalized_dict)

    panel = pd.Panel(dfs)
    mean_df = panel.mean(axis=0)
    sd_df = panel.std(axis=0)
    logger.debug("Mean shape : %s " % str(mean_df.shape))
    logger.debug("Sd  shape : %s " % str(sd_df.shape))

    return mean_df, sd_df

# ==============================================================
def calc_zscore_10mb_all(count_file, mean, sd):

    logger.debug("calc_zscore_10mb_all run")
    # -------------------------------------
    # 1. Calculate zscore by chromosome
    # -------------------------------------
    norm_dict = getNormalized10mbData(count_file, True)
    norm_XY_dict = getNormalized10mbData(count_file, False, key_start=len(norm_dict.keys()))
    logger.debug("norm_dict length %d" % len(norm_dict))
    logger.debug("norm_XY_dict length %d" % len(norm_XY_dict))

    merged_norm_dict = merge_two_dicts(norm_dict, norm_XY_dict)
    normalized_df = pd.DataFrame(merged_norm_dict)

    logger.debug("Merged Normalized DF shape : %s " % str(normalized_df.shape))

    z = {j : [(normalized_df.iloc[i,j] - mean.iloc[i,j])/sd.iloc[i,j] if sd.iloc[i,j] != 0.0 else 0.0 for i in range(len(mean))] for j in range(len(mean.iloc[0]))}

    zscore_10mb = pd.DataFrame(z)
    zscore_10mb.to_csv(count_file.split('10mb.')[0] + "zscore.10mb_all.txt", sep='\t')
    logger.debug(count_file.split('10mb.')[0] + "zscore.10mb_all.txt")

    return zscore_10mb

# ==============================================================
def draw_heatmap_10mb(zscore, ax):
    logger.debug("draw_heatmap_10mb run")
    x_ticks_loc = [0,  25,  50,  70,  90, 109, 127, 143, 158, 173, 187, 201, 215, 227, 238, 249, 259, 268, 276, 282, 289, 294, 300, 316]
    x_ticks_name = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    vline_loc = [x-0.5 for x in x_ticks_loc]
    #plt.vlines(vline_loc, -.5, 22)
    plt.xticks(x_ticks_loc, x_ticks_name, rotation=60)
    plt.yticks([x for x in range(22)],x_ticks_name);
    #ax.set_xticklabels(x_ticks_loc, rotation=70)
    #ax.set_yticklabels(x_ticks_name)
    ax.vlines(vline_loc, -.5, 22, colors='black', linestyles="solid", linewidth= 1.0)
    # Minor ticks
    ax.set_xticks(np.arange(-.5, 322, 1), minor=True);
    ax.set_yticks(np.arange(-.5, 22, 1), minor=True);
    
    ax.imshow(zscore, interpolation='none', cmap=plt.get_cmap('RdBu_r'), vmin=-7, vmax=7)
    ax.grid(which='minor', linestyle=':', linewidth=0.3)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

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
        logger.debug(filename_gendersearch)
        try:
            if gender_zipped_list[1][gender_zipped_list[0].index(filename_gendersearch)] == 'XX':
                zscore = calc_zscore_10mb_all(count_file, female_mean_df, female_sd_df)
                logger.debug("%s --> XX" %(filename_gendersearch))
            else:
                zscore = calc_zscore_10mb_all(count_file, male_mean_df, male_sd_df)
                logger.debug("%s --> XY" %(filename_gendersearch))

            if args.draw_fig == 'Y':
                ax_heatmap = plt.subplot(gs[2:,0:2])
                draw_heatmap_10mb(zscore, ax_heatmap)

                plt.tight_layout()
                filename_heatmap = os.path.basename(count_file).split(CFG_CountSplitForHeatmap)[0]
                plt.savefig("%s/%s.10mb_all.heatmap.png" %(CFG_SampleRunDir, filename_heatmap))
                logger.info("%s.10mb_all.heatmap.png created... OK" %(filename_heatmap))

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

    for zscore_file in zscore_list:

        ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
        zscore_10mb_file = '.'.join(zscore_file.split('.')[:-1]) + ".10mb.txt"
        zscore_10mb_df = pd.DataFrame.from_csv(zscore_10mb_file, sep='\t')
        draw_heatmap_10mb(zscore_10mb_df, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)

        plt.tight_layout()
        file_name = os.path.basename(zscore_file).split(CFG_CountSplitForHeatmap)[0]
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
        mean_df, sd_df = makeReferenceData(file_list, autosomal=True)
        #print mean_df.iloc[:,:5]

        # 2-2. make male reference
        logger.info("Generate male reference with %s" %(CFG_ReferenceMaleCountDir))
        file_list = glob.glob(CFG_ReferenceMaleCountDir+"/*"+CFG_CountExt)
        male_mean_df, male_sd_df = makeReferenceData(file_list, autosomal=False)
        #print male_mean_df[:,:5]

        # 2-3. make female reference
        logger.info("Generate female reference with %s" %(CFG_ReferenceFemaleCountDir))
        file_list = glob.glob(CFG_ReferenceFemaleCountDir+"/*"+CFG_CountExt)
        female_mean_df, female_sd_df = makeReferenceData(file_list, autosomal=False)
        #print female_mean_df[:,:5]

        male_mean_sum_df = pd.concat([mean_df, male_mean_df], axis=1)
        male_sd_sum_df = pd.concat([sd_df, male_sd_df], axis=1)
        female_mean_sum_df = pd.concat([mean_df, female_mean_df], axis=1)
        female_sd_sum_df = pd.concat([sd_df, female_sd_df], axis=1)

        # 2-4 overwrite csv files
        with open(CFG_MaleMean, 'w') as f:
            male_mean_sum_df.to_csv(f, sep='\t', index=False, header=False)
        with open(CFG_MaleSd, 'w') as f:
            male_sd_sum_df.to_csv(f, sep='\t', index=False, header=False)
        with open(CFG_FemaleMean, 'w') as f:
            female_mean_sum_df.to_csv(f, sep='\t', index=False, header=False)
        with open(CFG_FemaleSd, 'w') as f:
            female_sd_sum_df.to_csv(f, sep='\t', index=False, header=False)

        logger.debug(male_mean_sum_df.head())
        logger.debug("Male Mean DF : %s" % str(male_mean_sum_df.shape))
        logger.debug(female_sd_sum_df.head())
        logger.debug("Female Sd DF : %s" % str(female_sd_sum_df.shape))
        
    else:
        logger.info("Read reference files")
        male_mean_df = pd.read_csv(CFG_MaleMean, sep='\t', header=None)
        male_sd_df = pd.read_csv(CFG_MaleSd, sep='\t', header=None)
        female_mean_df = pd.read_csv(CFG_FemaleMean, sep='\t', header=None)
        female_sd_df = pd.read_csv(CFG_FemaleSd, sep='\t', header=None)

    # -------------------------------------------------------
    # 4. create normal & trisomy control z-score
    # -------------------------------------------------------
    plt.figure(dpi=200)
    plt.rcParams["figure.figsize"]= (28, 10)
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
                draw_heatmap(zscore, index_dict, color_list, ax_heatmap)

                logger.debug(zscore.head())

                ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
                draw_heatmap_10mb(zscore_10mb, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)
                logger.debug(zscore_10mb.head())

                plt.tight_layout()
                file_name = os.path.basename(count_file).split(CFG_CountSplitForHeatmap)[0]
                plt.savefig("%s/%s.heatmap.png" %(CFG_ReferenceMaleCountDir, file_name))
                logger.info("%s.heatmap.png created... OK" %(file_name))


        logger.info("Generate female normal z-score with %s" %(CFG_ReferenceFemaleCountDir))
        file_list = glob.glob(CFG_ReferenceFemaleCountDir+"/*"+CFG_CountExt)

        for count_file in file_list:
            zscore, index_dict, color_list, zscore_10mb, index_dict_10mb, color_list_10mb, row_count= calc_zscore_10mb(count_file, female_mean_df, female_sd_df, female_10mb_mean_df, female_10mb_sd_df)

            if args.draw_fig == 'Y':
                ax_heatmap = plt.subplot(gs[0:2,0:2])
                draw_heatmap(zscore, index_dict, color_list, ax_heatmap)

                ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
                draw_heatmap_10mb(zscore_10mb, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)

                plt.tight_layout()
                file_name = os.path.basename(count_file).split(CFG_CountSplitForHeatmap)[0]
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
                    draw_heatmap(zscore, index_dict, color_list, ax_heatmap)

                    logger.debug(zscore.head())

                    ax_heatmap_10mb = plt.subplot(gs[2:,0:2])
                    draw_heatmap_10mb(zscore_10mb, row_count, index_dict_10mb, color_list_10mb, ax_heatmap_10mb)
                    logger.debug(zscore_10mb.head())

                    plt.tight_layout()
                    file_name = os.path.basename(count_file).split(CFG_CountSplitForHeatmap)[0]
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
        logger.debug(filename_gendersearch)
        try:
            if gender_zipped_list[1][gender_zipped_list[0].index(filename_gendersearch)] == 'XX':
                zscore = calc_zscore_10mb_all(args.target_file, female_mean_df, female_sd_df)
                logger.debug("%s --> XX" %(filename_gendersearch))
            else:
                zscore = calc_zscore_10mb_all(args.target_file, male_mean_df, male_sd_df)
                logger.debug("%s --> XY" %(filename_gendersearch))

            if args.draw_fig == 'Y':
                ax_heatmap = plt.subplot(gs[2:,0:2])
                draw_heatmap_10mb(zscore, ax_heatmap)

                plt.tight_layout()
                #filename_heatmap = os.path.basename(args.target_file).split(CFG_CountSplitForHeatmap)[0]
                filename_heatmap = os.path.splitext(args.target_file)[0] + ".all.png"
                plt.savefig(filename_heatmap)
                logger.info("%s.10mb_all.heatmap.png created... OK" %(filename_heatmap))

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

