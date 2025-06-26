import pandas as pd

# hg38 염색체 크기 (단위: bp)
chr_sizes = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
    'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
    'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
    'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
    'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}

def calculate_coverage(bed_path):
    # BED 파일 로드
    df = pd.read_csv(bed_path, sep='\t', header=None, 
                     names=['chr', 'start', 'end'])
    
    # 각 영역 길이 계산
    df['length'] = df['end'] - df['start']
    
    # 염색체별 총 커버리지 계산
    chr_coverage = df.groupby('chr')['length'].sum().reset_index()
    
    # 전체 염색체 크기 정보 추가
    coverage_data = []
    for chr_name, total_size in chr_sizes.items():
        cov = chr_coverage[chr_coverage['chr'] == chr_name]['length'].sum()
        ratio = (cov / total_size * 100) if total_size > 0 else 0
        coverage_data.append({
            'Chromosome': chr_name,
            'Total Size (bp)': total_size,
            'Covered (bp)': cov,
            'Coverage (%)': round(ratio, 2)
        })
    
    return pd.DataFrame(coverage_data)

if __name__ == "__main__":
    df =calculate_coverage("/home/ken/ken-nipt/src/util/36mer/mappability_clean_bed/hg19_mappability_0.9_clean_all.bed")
    print(df)

