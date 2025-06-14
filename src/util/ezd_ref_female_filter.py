import pandas as pd
import sys
import os

def process_female_data_file(input_file):
    """
    Female 데이터 파일을 처리하는 함수
    """
    # 파일 읽기
    df = pd.read_csv(input_file, sep='\t')
    print(f"원본 데이터 행 수: {len(df)}")
    print(f"컬럼: {list(df.columns)}")
    
    # 1. XX/N/XXX 타입 중 ON2로 시작하는 샘플만 선택
    xx_n_xxx_pattern = r'^(XX|N|XXX)$'
    xx_n_xxx_mask = df['type'].str.contains(xx_n_xxx_pattern, na=False)
    on2_mask = df['Sample'].str.startswith('ON2', na=False)
    xx_filtered = df[xx_n_xxx_mask & on2_mask].copy()
    
    # N 타입을 XX로 변경
    xx_filtered.loc[xx_filtered['type'] == 'N', 'type'] = 'XX'
    
    # 2. XO 타입 모두 선택 (ON2 조건 없음)
    xo_pattern = r'^XO$'
    xo_mask = df['type'].str.contains(xo_pattern, na=False)
    xo_filtered = df[xo_mask].copy()
    
    # 3. TP/FP 타입 선택하고 괄호 제거
    tp_fp_pattern = r'\((TP|FP)\)'
    tp_fp_mask = df['type'].str.contains(tp_fp_pattern, na=False)
    tp_fp_filtered = df[tp_fp_mask].copy()
    tp_fp_filtered['type'] = tp_fp_filtered['type'].str.replace(r'\([^)]*\)', '', regex=True)
    
    print(f"XX/N/XXX 타입 중 ON2로 시작하는 데이터: {len(xx_filtered)}개")
    print(f"XO 타입 데이터: {len(xo_filtered)}개")
    print(f"TP/FP 타입 데이터: {len(tp_fp_filtered)}개")
    
    # 4. 모든 데이터 합치기
    final_df = pd.concat([xx_filtered, xo_filtered, tp_fp_filtered], ignore_index=True)
    
    # 5. Sample 컬럼 변경
    n_counter = 1
    p_counter = 1
    
    for idx in final_df.index:
        if final_df.loc[idx, 'type'] == 'XX':
            final_df.loc[idx, 'Sample'] = f'N{n_counter}'
            n_counter += 1
        else:
            final_df.loc[idx, 'Sample'] = f'P{p_counter}'
            p_counter += 1
    
    print(f"최종 데이터 행 수: {len(final_df)}")
    
    # 6. 파일 저장
    file_name, file_ext = os.path.splitext(input_file)
    output_file = f"{file_name}_filtered{file_ext}"
    final_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"결과 저장: {output_file}")
    
    # 결과 요약
    print("\n=== 처리 결과 ===")
    type_counts = final_df['type'].value_counts()
    for type_name, count in type_counts.items():
        print(f"{type_name} 타입: {count}개")
    
    print("\n샘플 미리보기:")
    print(final_df.head())
    
    return final_df

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("사용법: python script.py <입력파일경로>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if not os.path.exists(input_file):
        print(f"파일을 찾을 수 없습니다: {input_file}")
        sys.exit(1)
    
    try:
        result = process_female_data_file(input_file)
        print("처리 완료!")
    except Exception as e:
        print(f"오류 발생: {str(e)}")
        sys.exit(1)
