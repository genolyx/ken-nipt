import pyBigWig
import pybedtools
from tqdm import tqdm

# 입력 파일 경로
#mappability_bw = "wgEncodeCrgMapabilityAlign100mer.bw"
mappability_bw = "./36mer/wgEncodeCrgMapabilityAlign36mer.bw"
blacklist_bed = "ENCFF001TDO.bed"  # hg19 ENCODE blacklist BED
output_dir = "./36mer/mappability_clean_bed"

# 염색체 리스트
chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# 설정
step = 10000  # 10Kbp 단위 binning
score_cut = 1.0

# pyBigWig 열기
bw = pyBigWig.open(mappability_bw)
chrom_sizes = bw.chroms()

# 결과 저장
all_clean_bed = []

for chrom in tqdm(chroms, desc="Processing chromosomes"):
    if chrom not in chrom_sizes:
        continue

    chr_len = chrom_sizes[chrom]
    bins = []
    for start in range(0, chr_len, step):
        end = start + step
        if end > chr_len:
            end = chr_len
        try:
            score = bw.stats(chrom, start, end, type="mean")[0]
            if score is not None and score >= score_cut:
                bins.append([chrom, start, end])
        except Exception as e:
            print(f"Warning on {chrom}:{start}-{end} → {e}")
            continue

    # BED 파일로 저장
    raw_bed_file = f"./36mer/{chrom}_mappability_gt_{score_cut}.bed"
    with open(raw_bed_file, "w") as f:
        for b in bins:
            f.write(f"{b[0]}\t{b[1]}\t{b[2]}\n")

    # Blacklist 제외
    mappability_bed = pybedtools.BedTool(raw_bed_file)
    blacklist = pybedtools.BedTool(blacklist_bed)
    clean = mappability_bed.subtract(blacklist)

    # 결과 파일 저장
    clean_file = f"{output_dir}/{chrom}_mappability_clean.bed"
    clean.saveas(clean_file)

    all_clean_bed.append(clean)

bw.close()

# 전체 합치기
merged = pybedtools.BedTool.cat(*all_clean_bed, postmerge=False)
merged.saveas(f"{output_dir}/hg19_mappability_{score_cut}_clean_all.bed")

print("✅ All chromosomes processed and saved.")
