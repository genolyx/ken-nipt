#!/usr/bin/env python3
"""
PRIZM 10mb_all 파일만 재생성하는 스크립트
"""

import pandas as pd
import glob
import os
import sys
import re
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def natural_keys(text):
    def atoi(text):
        return int(text) if text.isdigit() else text
    return [atoi(c) for c in re.split(r'(\d+)', text)]

def getNormalized10mbData(filename, autosomal, key_start=0):
    df = pd.read_csv(filename, sep='\t')
    df = df[['chr', 'start', 'reads']].rename(columns={'start': 'bin', 'reads': 'count'})
    df['bin'] = df['bin'].apply(lambda x: (x-1)/10000000)
    df.bin = df.bin.astype(int)
    
    # Sort by chromosome order
    df = df.set_index('chr').loc[
        ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
         'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
         'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    ]
    df.reset_index(level=0, inplace=True)
    
    count_data_allsum = df.groupby(by=['chr'])['count'].sum()
    sum_df = count_data_allsum.to_frame()
    sorted_index = sorted(sum_df.index, key=natural_keys)
    sum_df = sum_df.loc[sorted_index]
    
    sum_df = sum_df.iloc[0:22, :]
    count_sum = sum_df['count'].sum()
    
    norm_data = [
        (sum_df.loc[sum_df.index[i], 'count']/float(count_sum))*100
        for i in range(len(sum_df))
    ]
    sum_df.insert(1, 'ratio', norm_data)
    
    if autosomal == True:
        df = df.loc[~df['chr'].isin(['chrX', 'chrY'])]
    else:
        df = df.loc[df['chr'].isin(['chrX', 'chrY'])]
    
    norm_10mb_data = [
        (df.loc[df.index[i], 'count']/float(count_sum))*100
        for i in range(len(df))
    ]
    df.insert(2, 'ratio', norm_10mb_data)
    
    normalized_dict = {
        j + key_start: [
            df.loc[df.index[j], 'ratio'] / sum_df.loc[sum_df.index[i], 'ratio']
            if sum_df.loc[sum_df.index[i], 'ratio'] != 0.0
            and df.loc[df.index[j], 'ratio'] != 0.0
            else 0
            for i in range(len(sum_df))
        ]
        for j in range(len(df))
    }
    
    return normalized_dict

def makeReference_10mball(file_list, autosomal):
    dfs = {}
    for idx, filename in enumerate(file_list):
        try:
            normalized_dict = getNormalized10mbData(filename, autosomal)
            dfs[idx] = pd.DataFrame.from_dict(normalized_dict)
        except Exception as e:
            logger.warning(f"Failed to process {filename}: {e}")
            continue
    
    if not dfs:
        return None, None
    
    all_dfs = pd.concat(dfs, axis=0, keys=range(len(dfs)))
    mean_df = all_dfs.groupby(level=1).mean()
    sd_df = all_dfs.groupby(level=1).std()
    
    return mean_df, sd_df

def collect_count_files(sample_df, group):
    """샘플 DataFrame으로부터 count 파일 수집"""
    count_files = []
    
    for _, row in sample_df.iterrows():
        sample_id = row['sample_id']
        sample_dir = row['sample_dir']
        
        # group에 따른 파일 이름 결정
        if group == 'orig':
            filename = f"{sample_id}.of_orig.10mb.wig.Normalization.txt"
        elif group == 'fetus':
            filename = f"{sample_id}.of_fetus.10mb.wig.Normalization.txt"
        elif group == 'mom':
            filename = f"{sample_id}.of_mom.10mb.wig.Normalization.txt"
        else:
            continue
        
        # Output_hmmcopy 디렉토리 안에 파일이 있음
        file_path = os.path.join(sample_dir, "Output_hmmcopy", filename)
        
        if os.path.exists(file_path):
            count_files.append(file_path)
    
    return count_files

def fix_10mb_all(group, sample_list_file, output_dir):
    """특정 group의 10mb_all 파일만 재생성"""
    
    logger.info(f"=== Fixing {group} 10mb_all files ===")
    
    # 샘플 리스트 읽기
    df = pd.read_csv(sample_list_file, sep='\t')
    
    # 성별 구분
    male_df = df[df['fetal_gender(gd_2)'] == 'XY']
    female_df = df[df['fetal_gender(gd_2)'] == 'XX']
    
    logger.info(f"Samples: Total={len(df)}, Male={len(male_df)}, Female={len(female_df)}")
    
    # Count 파일 수집
    logger.info(f"Collecting count files...")
    all_count_files = collect_count_files(df, group)
    male_count_files = collect_count_files(male_df, group)
    female_count_files = collect_count_files(female_df, group)
    
    logger.info(f"Count files found: Total={len(all_count_files)}, Male={len(male_count_files)}, Female={len(female_count_files)}")
    
    if len(all_count_files) == 0:
        logger.error(f"No count files found for {group}")
        return False
    
    # 1. Calculate autosomal bins (from all files)
    logger.info(f"Calculating autosomal 10mb_all...")
    total_mean_df, total_sd_df = makeReference_10mball(all_count_files, autosomal=True)
    
    if total_mean_df is None:
        logger.error("Failed to calculate autosomal 10mb_all reference")
        return False
    
    logger.info(f"  Autosomal shape: {total_mean_df.shape}")
    
    # 2. Calculate male sex chromosome bins
    logger.info(f"Calculating male sex chromosome 10mb_all...")
    male_sex_mean_df, male_sex_sd_df = makeReference_10mball(male_count_files, autosomal=False)
    
    if male_sex_mean_df is not None:
        logger.info(f"  Male sex chr shape: {male_sex_mean_df.shape}")
        
        # Combine autosomal + sex chromosomes
        male_mean_10mball_df = pd.concat([total_mean_df, male_sex_mean_df], axis=1)
        male_sd_10mball_df = pd.concat([total_sd_df, male_sex_sd_df], axis=1)
        
        male_mean_file = os.path.join(output_dir, 'male_10mb_all_mean.csv')
        male_sd_file = os.path.join(output_dir, 'male_10mb_all_sd.csv')
        
        male_mean_10mball_df.to_csv(male_mean_file, sep='\t', index=False, header=False)
        male_sd_10mball_df.to_csv(male_sd_file, sep='\t', index=False, header=False)
        
        logger.info(f"  ✓ Saved: {male_mean_file} (shape: {male_mean_10mball_df.shape})")
        logger.info(f"  ✓ Saved: {male_sd_file}")
    
    # 3. Calculate female sex chromosome bins
    logger.info(f"Calculating female sex chromosome 10mb_all...")
    female_sex_mean_df, female_sex_sd_df = makeReference_10mball(female_count_files, autosomal=False)
    
    if female_sex_mean_df is not None:
        logger.info(f"  Female sex chr shape: {female_sex_mean_df.shape}")
        
        # Combine autosomal + sex chromosomes
        female_mean_10mball_df = pd.concat([total_mean_df, female_sex_mean_df], axis=1)
        female_sd_10mball_df = pd.concat([total_sd_df, female_sex_sd_df], axis=1)
        
        female_mean_file = os.path.join(output_dir, 'female_10mb_all_mean.csv')
        female_sd_file = os.path.join(output_dir, 'female_10mb_all_sd.csv')
        
        female_mean_10mball_df.to_csv(female_mean_file, sep='\t', index=False, header=False)
        female_sd_10mball_df.to_csv(female_sd_file, sep='\t', index=False, header=False)
        
        logger.info(f"  ✓ Saved: {female_mean_file} (shape: {female_mean_10mball_df.shape})")
        logger.info(f"  ✓ Saved: {female_sd_file}")
    
    return True

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Fix PRIZM 10mb_all files')
    parser.add_argument('--sample-list', required=True, help='Filtered sample list TSV')
    parser.add_argument('--group', required=True, choices=['orig', 'fetus', 'mom'], help='Group to fix')
    parser.add_argument('--output-dir', required=True, help='Output directory (e.g., /path/to/ucl_new/PRIZM/orig)')
    
    args = parser.parse_args()
    
    success = fix_10mb_all(args.group, args.sample_list, args.output_dir)
    
    if success:
        logger.info("✅ Successfully fixed 10mb_all files!")
        sys.exit(0)
    else:
        logger.error("❌ Failed to fix 10mb_all files")
        sys.exit(1)
