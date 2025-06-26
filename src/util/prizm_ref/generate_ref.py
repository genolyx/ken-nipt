#!/usr/bin/miniconda3/bin/python3
"""
---------------------------------------------
Generate PRIZM reference

Author: {author}
Contact: {email}
---------------------------------------------
"""

import glob
import pandas as pd
import numpy as np
import sys
import os
import json
import logging
import re
from pathlib import Path

__author__ = 'Kenneth Kwon'
__email__ = "kennethkwon@genecurate.xyz"
__version__ = '0.9'

# ==============================================================
class Config:

    REQUIRED_OUTPUT_FILES = [
        'TotalMean', 'TotalSd', 'Total10mbMean', 'Total10mbSd',
        'MaleMean', 'MaleSd', 'Male10mbMean', 'Male10mbSd',
        'FemaleMean', 'FemaleSd', 'Female10mbMean', 'Female10mbSd'
    ]

    def __init__(self, config_path, reference_type='orig'):
        with open(config_path, 'r') as f:
            config = json.load(f)

        self._validate_config(config, reference_type)

        # Set base directories
        self.input_dir = config['paths']['input_dirs'][reference_type]
        self.output_dir = config['paths']['output_dirs'][reference_type]

        # Set gender-specific input directories
        self.male_input_dir = os.path.join(self.input_dir, 'M')
        self.female_input_dir = os.path.join(self.input_dir, 'F')

        # Set output file paths
        output_files = config['paths']['output_files']
        for file_key, file_name in output_files.items():
            full_path = os.path.join(self.output_dir, file_name)
            setattr(self, file_key, full_path)

        self._create_output_dirs()

    def _validate_config(self, config, reference_type):
        """Validate that all required paths are in the config"""
        if 'paths' not in config:
            raise ValueError("Configuration must contain 'paths' section")

        # Validate input_dirs and output_dirs
        for dir_type in ['input_dirs', 'output_dirs']:
            if dir_type not in config['paths']:
                raise ValueError(f"Configuration must contain '{dir_type}' section")
            if reference_type not in config['paths'][dir_type]:
                raise ValueError(f"Missing {reference_type} in {dir_type}")

        # Validate output_files
        if 'output_files' not in config['paths']:
            raise ValueError("Configuration must contain 'output_files' section")

        for file in self.REQUIRED_OUTPUT_FILES:
            if file not in config['paths']['output_files']:
                raise ValueError(f"Missing required output file: {file}")

    def _create_output_dirs(self):
        """Create output directory if it doesn't exist"""
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

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
    return [atoi(c) for c in re.split(r'(\d+)', text)]

# ==============================================================
def getNormalized10mbData(filename, autosomal, key_start=0):

    df = pd.read_csv(filename, sep='\t')
    df = df[['chr', 'start', 'reads']].rename(columns={'start':'bin', 'reads':'count'})
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

    normalized_dict = {j+key_start : [df.loc[df.index[j],'ratio']/sum_df.loc[sum_df.index[i],'ratio'] if sum_df.loc[sum_df.index[i], 'ratio'] != 0.0 and df.loc[df.index[j], 'ratio'] != 0.0 else 0 for i in range(len(sum_df))] for j in range(len(df))}

    return normalized_dict

# ==============================================================
def makeReference_10mball(file_list, autosomal):

    dfs = {}
    for idx, filename in enumerate(file_list):
        normalized_dict = getNormalized10mbData(filename, autosomal)
        dfs[idx] = pd.DataFrame.from_dict(normalized_dict)

    # Stack all DataFrames into a single 3D structure
    all_dfs = pd.concat(dfs, axis=0, keys=range(len(dfs)))

    # Calculate mean and standard deviation across the first axis
    mean_df = all_dfs.groupby(level=1).mean()
    sd_df = all_dfs.groupby(level=1).std()

    logger.debug("Mean shape : %s " % str(mean_df.shape))
    logger.debug("Sd  shape : %s " % str(sd_df.shape))

    return mean_df, sd_df

# ==============================================================
def makeReference_10mb(file_list, mean_file, sd_file, mean_10mb_file, sd_10mb_file, autosomal):
    # 1. Reference by chromosome count
    dfs = {}
    for idx, filename in enumerate(file_list):
        count_data = pd.read_csv(filename, sep='\t')
        count_data = count_data[['chr', 'start', 'reads']].rename(columns={'start':'bin', 'reads':'count'})

        count_data['bin'] = count_data['bin'].apply(lambda x: x/10000000)
        count_data_allsum = count_data.groupby(by=['chr'])['count'].sum()
        count_data_df = count_data_allsum.to_frame()
        sorted_index = sorted(count_data_df.index, key=natural_keys)
        count_data_df = count_data_df.loc[sorted_index]

        if autosomal:
            count_data_df = count_data_df.iloc[0:22, :]

        count_sum = sum(count_data_df.loc[:,'count'])
        norm_data = [(count_data_df.iloc[i,count_data_df.columns.get_loc('count')]/float(count_sum))*100 for i in range(len(count_data_df))]

        count_data_df.insert(1, 'ratio', norm_data)
        normalized_dict = {j : [count_data_df.iloc[i,count_data_df.columns.get_loc('ratio')]/count_data_df.iloc[j,count_data_df.columns.get_loc('ratio')] for i in range(len(count_data_df))] for j in range(len(count_data_df))}

        #normalized_df = pd.DataFrame(normalized_dict)

        dfs[idx] = pd.DataFrame.from_dict(normalized_dict)

    # Replace Panel with concat and stack
    all_dfs = pd.concat(dfs, axis=0)
    mean_df = all_dfs.groupby(level=1).mean()
    sd_df = all_dfs.groupby(level=1).std()

    mean_df.to_csv(mean_file, sep='\t', index=False, header=False)
    sd_df.to_csv(sd_file, sep='\t', index=False, header=False)

    # 2. Reference by 10mb bin count
    dfs = {}
    for idx, filename in enumerate(file_list):
        count_data = pd.read_csv(filename, sep='\t')
        count_data = count_data[['chr', 'start', 'reads']].rename(columns={'start':'bin', 'reads':'count'})

        count_data['bin'] = count_data['bin'].apply(lambda x: x/10000000)
        count_data_allsum = count_data.groupby(by=['chr'])['count'].sum()
        count_data_df = count_data_allsum.to_frame()
        count_data_df.insert(0, 'bin', 0)
        sorted_index = sorted(count_data_df.index, key=natural_keys)

        count_data_df = count_data_df.drop(['chr9', 'chr13', 'chr16', 'chr18', 'chr21', 'chr22', 'chrX', 'chrY'])
        excluded_chromosomes = ['chr9', 'chr13', 'chr16', 'chr18', 'chr21', 'chr22']
        if not autosomal:
            excluded_chromosomes.extend(['chrX', 'chrY'])
            
        count_data_10mb = count_data.set_index('chr').loc[excluded_chromosomes]

        count_data_mix = pd.concat([count_data_df, count_data_10mb])
        count_data_mix = count_data_mix.loc[sorted_index[:-2] if autosomal else sorted_index]

        count_sum = count_data_mix.loc[:,'count'].sum()
        norm_data = [(count_data_mix.iloc[i,count_data_mix.columns.get_loc('count')]/float(count_sum))*100 for i in range(len(count_data_mix))]
        count_data_mix.insert(2, 'ratio', norm_data)

        normalized_dict = {j : [count_data_mix.iloc[i,count_data_mix.columns.get_loc('ratio')]/count_data_mix.iloc[j,count_data_mix.columns.get_loc('ratio')] if count_data_mix.iloc[i, count_data_mix.columns.get_loc('ratio')] != 0.0 and count_data_mix.iloc[j, count_data_mix.columns.get_loc('ratio')] != 0.0 else 0 for i in range(len(count_data_mix))] for j in range(len(count_data_mix))}

        #normalized_df = pd.DataFrame(normalized_dict)

        dfs[idx] = pd.DataFrame.from_dict(normalized_dict)

    # Replace Panel with concat and stack
    all_10mb_dfs = pd.concat(dfs, axis=0)
    mean_10mb_df = all_10mb_dfs.groupby(level=1).mean()
    sd_10mb_df = all_10mb_dfs.groupby(level=1).std()

    mean_10mb_df.to_csv(mean_10mb_file, sep='\t', index=False, header=False)
    sd_10mb_df.to_csv(sd_10mb_file, sep='\t', index=False, header=False)

    return mean_df, sd_df, mean_10mb_df, sd_10mb_df

def generate_references(cfg, postfix):
    """
    Generate reference files from male and female data
    """
    # 2-1. Make overall reference using both M and F directories
    logger.info(cfg.male_input_dir)
    logger.info(cfg.female_input_dir)

    if postfix == "*10mb.txt":
        male_files = [f for f in glob.glob(os.path.join(cfg.male_input_dir, postfix))
                    if not ('filter' in f or 'filter1' in f or 'filter_out' in f)]
        female_files = [f for f in glob.glob(os.path.join(cfg.female_input_dir, postfix))
                    if not ('filter' in f or 'filter1' in f or 'filter_out' in f)]
    else:
        # For other patterns, use exact match
        male_files = glob.glob(os.path.join(cfg.male_input_dir, postfix))
        female_files = glob.glob(os.path.join(cfg.female_input_dir, postfix))

    all_files = male_files + female_files

    print(f"Generating total reference from {len(all_files)} files")
    total_mean_df, total_sd_df, total_10mb_mean_df, total_10mb_sd_df = makeReference_10mb(
        all_files,
        cfg.TotalMean,
        cfg.TotalSd,
        cfg.Total10mbMean,
        cfg.Total10mbSd,
        autosomal=True
    )

    # 2-2. Generate male reference
    print(f"Generating male reference from {len(male_files)} files")
    male_mean_df, male_sd_df, male_10mb_mean_df, male_10mb_sd_df = makeReference_10mb(
        male_files,
        cfg.MaleMean,
        cfg.MaleSd,
        cfg.Male10mbMean,
        cfg.Male10mbSd,
        autosomal=False
    )

    # Generate female reference
    print(f"Generating female reference from {len(female_files)} files")
    female_mean_df, female_sd_df, female_10mb_mean_df, female_10mb_sd_df = makeReference_10mb(
        female_files,
        cfg.FemaleMean,
        cfg.FemaleSd,
        cfg.Female10mbMean,
        cfg.Female10mbSd,
        autosomal=False
    )

    # 2-3. Update gender-specific references with total reference data
    for gender_mean_df in [male_mean_df, female_mean_df]:
        gender_mean_df.update(total_mean_df)

    for gender_sd_df in [male_sd_df, female_sd_df]:
        gender_sd_df.update(total_sd_df)

    for gender_10mb_mean_df in [male_10mb_mean_df, female_10mb_mean_df]:
        gender_10mb_mean_df.update(total_10mb_mean_df)

    for gender_10mb_sd_df in [male_10mb_sd_df, female_10mb_sd_df]:
        gender_10mb_sd_df.update(total_10mb_sd_df)

    # 2-4. Save updated references
    output_pairs = [
        (male_mean_df, cfg.MaleMean),
        (male_sd_df, cfg.MaleSd),
        (female_mean_df, cfg.FemaleMean),
        (female_sd_df, cfg.FemaleSd),
        (male_10mb_mean_df, cfg.Male10mbMean),
        (male_10mb_sd_df, cfg.Male10mbSd),
        (female_10mb_mean_df, cfg.Female10mbMean),
        (female_10mb_sd_df, cfg.Female10mbSd)
    ]

    for df, filepath in output_pairs:
        print(f"Saving to {filepath}")
        with open(filepath, 'w') as f:
            df.to_csv(f, sep='\t', index=False, header=False)

def generate_10mball_references(cfg, postfix):

    # 2-1. Make overall reference using both M and F directories
    logger.info(cfg.male_input_dir)
    logger.info(cfg.female_input_dir)

    if postfix == "*10mb.txt":
        male_files = [f for f in glob.glob(os.path.join(cfg.male_input_dir, postfix))
                    if not ('filter' in f or 'filter1' in f or 'filter_out' in f)]
        female_files = [f for f in glob.glob(os.path.join(cfg.female_input_dir, postfix))
                    if not ('filter' in f or 'filter1' in f or 'filter_out' in f)]
    else:
        # For other patterns, use exact match
        male_files = glob.glob(os.path.join(cfg.male_input_dir, postfix))
        female_files = glob.glob(os.path.join(cfg.female_input_dir, postfix))

    all_files = male_files + female_files

    print(f"Generating total reference from {len(all_files)} files")
    mean_df, sd_df = makeReference_10mball(all_files, autosomal=True)

    # 2-2. Generate male reference
    print(f"Generating male reference from {len(male_files)} files")
    male_mean_df, male_sd_df = makeReference_10mball(male_files, autosomal=False)

    # Generate female reference
    print(f"Generating female reference from {len(female_files)} files")
    female_mean_df, female_sd_df = makeReference_10mball(female_files, autosomal=False)

    male_mean_10mball_df = pd.concat([mean_df, male_mean_df], axis=1)
    male_sd_10mball_df = pd.concat([sd_df, male_sd_df], axis=1)
    female_mean_10mball_df = pd.concat([mean_df, female_mean_df], axis=1)
    female_sd_10mball_df = pd.concat([sd_df, female_sd_df], axis=1)

    # 2-4. Save updated references
    output_pairs = [
        (male_mean_10mball_df, cfg.Male10mbAllMean),
        (male_sd_10mball_df, cfg.Male10mbAllSd),
        (female_mean_10mball_df, cfg.Female10mbAllMean),
        (female_sd_10mball_df, cfg.Female10mbAllSd),
    ]

    for df, filepath in output_pairs:
        print(f"Saving to {filepath}")
        with open(filepath, 'w') as f:
            df.to_csv(f, sep='\t', index=False, header=False)


# ==============================================================
#   MAIN FUNCTION
# ==============================================================
if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    config_path = "/Work/NIPT/refs/cordlife_new/PRIZM/PRIZM.json"

    try:
        # Create config objects for different reference types
        orig_cfg = Config(config_path, reference_type='orig')
        fetus_cfg = Config(config_path, reference_type='fetus')
        fetus1_cfg = Config(config_path, reference_type='fetus1')
        mom_cfg = Config(config_path, reference_type='mom')

        logger.info(f"Generating references from {orig_cfg.input_dir}")
        generate_references(orig_cfg, "*10mb.txt")
        generate_10mball_references(orig_cfg, "*10mb.txt")

        logger.info(f"Generating references from {fetus_cfg.input_dir}")
        generate_references(fetus_cfg, "*.filter.10mb.txt")
        generate_10mball_references(fetus_cfg, "*.filter.10mb.txt")

        #logger.info(f"Generating references from {fetus1_cfg.input_dir}")
        generate_references(fetus1_cfg, "*.filter1.10mb.txt")

        logger.info(f"Generating references from {mom_cfg.input_dir}")
        generate_references(mom_cfg, "*filter_out.10mb.txt")
        generate_10mball_references(mom_cfg, "*.filter_out.10mb.txt")

    except Exception as e:
        logger.error(f"Error: {str(e)}")
        raise

