import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import sys
import os

def load_male_data(file_path):
    """
    male.txt 파일 로드
    
    Args:
        file_path (str): male.txt 파일 경로
    
    Returns:
        pandas.DataFrame: 로드된 데이터
    """
    try:
        # 컬럼명이 UAR.X, UAR.Y 형태이므로 적절히 처리
        df = pd.read_csv(file_path, sep='\t')
        
        # 컬럼명 정리 (점을 언더스코어로 변경)
        df.columns = [col.replace('.', '_') for col in df.columns]
        
        print(f"데이터 로드 완료: {len(df)}개 샘플")
        print(f"컬럼: {list(df.columns)}")
        print(f"타입별 분포:")
        print(df['type'].value_counts())
        
        return df
    
    except Exception as e:
        print(f"파일 로드 오류: {e}")
        return None

def extract_xy_xyy_data(df):
    """
    XY(정상)와 XYY 타입 데이터 추출
    
    Args:
        df (pandas.DataFrame): 전체 데이터
    
    Returns:
        tuple: (xy_data, xyy_data)
    """
    xy_data = df[df['type'] == 'XY'].copy()
    xyy_data = df[df['type'] == 'XYY'].copy()
    
    print(f"\nXY (정상) 샘플: {len(xy_data)}개")
    print(f"XYY 샘플: {len(xyy_data)}개")
    
    if len(xy_data) == 0:
        print("Warning: XY 타입 데이터가 없습니다!")
    if len(xyy_data) == 0:
        print("Warning: XYY 타입 데이터가 없습니다!")
    
    return xy_data, xyy_data

def calculate_boundary_line(xy_data, percentile=95, ur_x_col='UAR_X', ur_y_col='UAR_Y'):
    """
    XY 데이터에서 경계선 계산
    
    Args:
        xy_data (pandas.DataFrame): XY 타입 데이터
        percentile (int): 경계선으로 사용할 백분위수 (기본: 95)
        ur_x_col (str): UAR_X 컬럼명
        ur_y_col (str): UAR_Y 컬럼명
    
    Returns:
        tuple: (slope, intercept, boundary_points)
    """
    if len(xy_data) == 0:
        print("XY 데이터가 없어서 경계선을 계산할 수 없습니다.")
        return None, None, []
    
    # UAR_X 범위 확인
    x_min, x_max = xy_data[ur_x_col].min(), xy_data[ur_x_col].max()
    print(f"\nXY 데이터 UAR_X 범위: {x_min:.3f} ~ {x_max:.3f}")
    print(f"XY 데이터 UAR_Y 범위: {xy_data[ur_y_col].min():.3f} ~ {xy_data[ur_y_col].max():.3f}")
    
    # UAR_X를 구간별로 나누어 각 구간의 상위 백분위수 구하기
    n_bins = min(10, len(xy_data) // 5)  # 구간 수는 데이터 양에 따라 조정
    x_bins = np.linspace(x_min, x_max, n_bins + 1)
    
    boundary_points = []
    
    for i in range(len(x_bins) - 1):
        x_bin_min, x_bin_max = x_bins[i], x_bins[i + 1]
        x_center = (x_bin_min + x_bin_max) / 2
        
        # 해당 구간의 데이터
        bin_mask = ((xy_data[ur_x_col] >= x_bin_min) & 
                   (xy_data[ur_x_col] < x_bin_max))
        bin_data = xy_data[bin_mask]
        
        if len(bin_data) >= 3:  # 최소 3개 이상의 데이터가 있어야 신뢰성 있음
            y_boundary = np.percentile(bin_data[ur_y_col], percentile)
            boundary_points.append((x_center, y_boundary))
            print(f"구간 [{x_bin_min:.2f}, {x_bin_max:.2f}]: {len(bin_data)}개 샘플, {percentile}th percentile = {y_boundary:.4f}")
    
    if len(boundary_points) < 2:
        print("경계선 계산을 위한 충분한 데이터가 없습니다.")
        return None, None, boundary_points
    
    # 선형 회귀로 경계선 구하기
    X = np.array([p[0] for p in boundary_points]).reshape(-1, 1)
    y = np.array([p[1] for p in boundary_points])
    
    reg = LinearRegression().fit(X, y)
    slope = reg.coef_[0]
    intercept = reg.intercept_
    
    # R² 점수 계산
    r2_score = reg.score(X, y)
    
    print(f"\n=== 경계선 계산 결과 ===")
    print(f"경계선 방정식: y = {slope:.6f}x + {intercept:.6f}")
    print(f"R² 점수: {r2_score:.4f}")
    print(f"사용된 구간 수: {len(boundary_points)}개")
    
    return slope, intercept, boundary_points

def test_xyy_detection(slope, intercept, test_x, test_y):
    """
    특정 좌표에서 XYY 여부 판정
    
    Args:
        slope (float): 경계선 기울기
        intercept (float): 경계선 y절편
        test_x (float): 테스트 UAR_X 값
        test_y (float): 테스트 UAR_Y 값
    
    Returns:
        str: 판정 결과
    """
    if slope is None or intercept is None:
        return "경계선이 계산되지 않음"
    
    boundary_y = slope * test_x + intercept
    margin = 0.005  # 여유 마진
    
    if test_y > boundary_y + margin:
        return "XYY Detected"
    elif test_y > boundary_y:
        return "XYY Suspected"
    else:
        return "Normal"

def plot_boundary_analysis(df, slope, intercept, boundary_points, ur_x_col='UAR_X', ur_y_col='UAR_Y'):
    """
    경계선 분석 결과 시각화
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # XY (정상) 데이터
    xy_data = df[df['type'] == 'XY']
    if len(xy_data) > 0:
        ax.scatter(xy_data[ur_x_col], xy_data[ur_y_col], 
                  color='blue', alpha=0.6, s=20, label=f'XY (Normal) - {len(xy_data)}개')
    
    # XYY 데이터
    xyy_data = df[df['type'] == 'XYY']
    if len(xyy_data) > 0:
        ax.scatter(xyy_data[ur_x_col], xyy_data[ur_y_col], 
                  color='green', alpha=0.8, s=40, marker='^', label=f'XYY - {len(xyy_data)}개')
    
    # 기타 타입들
    other_types = df[~df['type'].isin(['XY', 'XYY'])]
    for type_name in other_types['type'].unique():
        type_data = other_types[other_types['type'] == type_name]
        ax.scatter(type_data[ur_x_col], type_data[ur_y_col], 
                  alpha=0.7, s=30, label=f'{type_name} - {len(type_data)}개')
    
    # 경계선 구간별 포인트
    if boundary_points:
        bp_x = [p[0] for p in boundary_points]
        bp_y = [p[1] for p in boundary_points]
        ax.scatter(bp_x, bp_y, color='red', s=80, marker='x', 
                  linewidth=3, label='95th Percentile Points')
    
    # 경계선 그리기
    if slope is not None and intercept is not None:
        x_range = np.linspace(df[ur_x_col].min() * 0.95, df[ur_x_col].max() * 1.05, 100)
        y_boundary = slope * x_range + intercept
        ax.plot(x_range, y_boundary, 
               color='red', linestyle='-', linewidth=3, 
               label=f'XYY Boundary: y={slope:.4f}x+{intercept:.4f}')
        
        # XYY 영역 표시
        ax.fill_between(x_range, y_boundary, ax.get_ylim()[1], 
                       color='red', alpha=0.1, label='XYY Detection Zone')
    
    ax.set_xlabel('UAR[X] (%)', fontsize=12)
    ax.set_ylabel('UAR[Y] (%)', fontsize=12)
    ax.set_title('XY-XYY Boundary Analysis', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig, ax

def main(male_file_path):
    """
    메인 함수
    
    Args:
        male_file_path (str): male.txt 파일 경로
    """
    print("=== Male 데이터 XYY 경계선 분석 ===")
    
    # 1. 데이터 로드
    df = load_male_data(male_file_path)
    if df is None:
        return
    
    # 2. XY, XYY 데이터 추출
    xy_data, xyy_data = extract_xy_xyy_data(df)
    
    # 3. 경계선 계산
    slope, intercept, boundary_points = calculate_boundary_line(xy_data)
    
    # 4. 결과 출력
    if slope is not None:
        print(f"\n=== 최종 결과 ===")
        print(f"추천 경계선: y = {slope:.6f}x + {intercept:.6f}")
        print("이 선 위쪽에 있으면 XYY로 판정됩니다.")
        
        # 5. XYY 샘플들로 검증
        if len(xyy_data) > 0:
            print(f"\n=== XYY 샘플 검증 ===")
            correct_count = 0
            for idx, row in xyy_data.iterrows():
                result = test_xyy_detection(slope, intercept, row['UAR_X'], row['UAR_Y'])
                is_correct = "XYY" in result
                correct_count += is_correct
                print(f"{row['Sample']}: UAR_X={row['UAR_X']:.3f}, UAR_Y={row['UAR_Y']:.3f} → {result} {'✓' if is_correct else '✗'}")
            
            accuracy = correct_count / len(xyy_data) * 100
            print(f"XYY 검출 정확도: {accuracy:.1f}% ({correct_count}/{len(xyy_data)})")
    
    # 6. 시각화
    fig, ax = plot_boundary_analysis(df, slope, intercept, boundary_points)
    plt.show()
    
    return slope, intercept

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("사용법: python script.py <male.txt 파일 경로>")
        print("예시: python script.py /path/to/male.txt")
        sys.exit(1)
    
    male_file = sys.argv[1]
    
    if not os.path.exists(male_file):
        print(f"파일을 찾을 수 없습니다: {male_file}")
        sys.exit(1)
    
    slope, intercept = main(male_file)
    
    if slope is not None:
        print(f"\n코드에서 사용할 값:")
        print(f"slope = {slope}")
        print(f"intercept = {intercept}")
