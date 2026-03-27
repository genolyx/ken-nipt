import pandas as pd
import numpy as np

df = pd.read_csv("/home/ken/ken-nipt/reference_sample_list_GNCI.tsv", sep="\t")

xy = df[df["fetal_gender(gd_2)"]=="XY"][["YFF_2","SeqFF"]].dropna()

# IQR로 아주 튀는 값만 제거(권장)
def iqr_mask(s, k=1.5):
    q1, q3 = np.percentile(s, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - k*iqr, q3 + k*iqr
    return (s >= lo) & (s <= hi)

mask = iqr_mask(xy["YFF_2"]) & iqr_mask(xy["SeqFF"])
xyf = xy[mask]

# 선형회귀: YFF_2 = a + b*SeqFF
b, a = np.polyfit(xyf["SeqFF"].values, xyf["YFF_2"].values, 1)
print("a(intercept)=", a, " b(slope)=", b)

# 보정 적용
df["SeqFF_corrected"] = a + b*df["SeqFF"]
df["SeqFF_corrected"] = df["SeqFF_corrected"].clip(lower=0, upper=100)

df.to_csv("/home/ken/ken-nipt/reference_sample_list.with_seqff_corrected.tsv", sep="\t", index=False)
