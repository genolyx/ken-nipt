import pandas as pd
import sys
import os

def process_data_file(input_file):
    """
    데이터 파일을 처리하는 함수
    
    Args:
        input_file (str): 입력 파일 경로
    
    Returns:
        pandas.DataFrame: 처리된 데이터프레임
    """
    try:
        # 파일 읽기 (탭으로 구분된 것으로 가정)
        df = pd.read_csv(input_file, sep='\t')
        
        print(f"원본 데이터 행 수: {len(df)}")
        print(f"컬럼: {list(df.columns)}")
        
        # 1. type 컬럼에서 10P, TP, FP가 포함된 것만 필터링
        pattern = r'\((10P|TP|FP)\)'
        filtered_df = df[df['type'].str.contains(pattern, na=False)].copy()
        
        print(f"필터링 후 데이터 행 수: {len(filtered_df)}")
        
        # 2. Sample 컬럼을 P1, P2, P3... 형태로 변경
        for i, idx in enumerate(filtered_df.index, 1):
            filtered_df.loc[idx, 'Sample'] = f'P{i}'
        
        # 3. type 컬럼에서 괄호와 괄호 안의 내용 제거 (XXY(10P) -> XXY)
        filtered_df['type'] = filtered_df['type'].str.replace(r'\([^)]*\)', '', regex=True)
        
        # 4. 원본 데이터에서 type이 XY인 것들 처리
        xy_pattern = r'^XY$'  # 정확히 XY인 것만
        xy_mask = df['type'].str.contains(xy_pattern, na=False)
        
        if xy_mask.any():
            # XY 타입 중에서 Sample이 "ON2"로 시작하는 것만 선택
            xy_df = df[xy_mask & df['Sample'].str.startswith('ON2', na=False)].copy()
            
            if len(xy_df) > 0:
                # Sample 컬럼을 N1, N2, N3... 형태로 변경
                for i, idx in enumerate(xy_df.index, 1):
                    xy_df.loc[idx, 'Sample'] = f'N{i}'
                # type은 XY 그대로 유지
                
                print(f"XY 타입 중 ON2로 시작하는 데이터 행 수: {len(xy_df)}")
                
                # 필터링된 데이터와 XY 데이터를 합치기
                final_df = pd.concat([filtered_df, xy_df], ignore_index=True)
            else:
                print("XY 타입 중 ON2로 시작하는 샘플이 없습니다.")
                final_df = filtered_df
        else:
            print("XY 타입 데이터를 찾을 수 없습니다.")
            final_df = filtered_df
        
        print(f"최종 데이터 행 수: {len(final_df)}")
        
        # 5. 출력 파일명 생성 (_filtered 추가)
        file_name, file_ext = os.path.splitext(input_file)
        output_file = f"{file_name}_filtered{file_ext}"
        
        # 파일 저장
        final_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"결과가 저장되었습니다: {output_file}")
        
        # 결과 요약 출력
        print("\n=== 처리 결과 요약 ===")
        type_counts = final_df['type'].value_counts()
        for type_name, count in type_counts.items():
            print(f"{type_name} 타입 샘플 수: {count}")
        
        print("\n샘플 미리보기:")
        print(final_df.head(10))
        
        return final_df
        
    except Exception as e:
        print(f"오류 발생: {str(e)}")
        return None

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("사용법: python script.py <입력파일경로>")
        print("예시: python script.py data.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if not os.path.exists(input_file):
        print(f"파일을 찾을 수 없습니다: {input_file}")
        sys.exit(1)
    
    result = process_data_file(input_file)
    
    if result is not None:
        print("처리가 완료되었습니다!")
    else:
        print("처리 중 오류가 발생했습니다.")
        sys.exit(1)
