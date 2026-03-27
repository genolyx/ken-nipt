#!/usr/bin/env python3
"""
Case-insensitive 필터링 테스트 스크립트

다양한 대소문자 조합의 Result/MDResult 값이 제대로 필터링되는지 테스트
"""

import pandas as pd
import sys
import os

# Add the script directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_case_insensitive_filtering():
    """Case-insensitive 필터링 테스트"""
    
    print("="*60)
    print("Case-Insensitive Filtering Test")
    print("="*60)
    
    # 테스트 데이터 생성 (다양한 대소문자 조합 + 오매칭 방지)
    test_data = {
        'sample_id': [
            'SAMPLE001', 'SAMPLE002', 'SAMPLE003', 'SAMPLE004',
            'SAMPLE005', 'SAMPLE006', 'SAMPLE007', 'SAMPLE008',
            'SAMPLE009', 'SAMPLE010', 'SAMPLE011', 'SAMPLE012',
            'SAMPLE013', 'SAMPLE014', 'SAMPLE015'
        ],
        'Result': [
            'Low Risk',      # OK
            'High Risk',     # 제외 (정확한 표기)
            'high risk',     # 제외 (소문자)
            'High risk',     # 제외 (혼합)
            'LOW RISK',      # OK
            'No Call',       # 제외 (정확한 표기)
            'no call',       # 제외 (소문자)
            'No call',       # 제외 (혼합)
            'Low Risk',      # OK
            'HIGH RISK',     # 제외 (대문자)
            'NO CALL',       # 제외 (대문자)
            'Low Risk',      # OK
            'None',          # OK (부분 매칭 방지 테스트)
            'Recalled',      # OK (부분 매칭 방지 테스트)
            'Normal'         # OK (부분 매칭 방지 테스트)
        ],
        'MDResult': [
            'Low Risk',      # OK
            'Low Risk',      # OK
            'High Risk',     # 제외
            'No Call',       # 제외
            'Low Risk',      # OK
            'Low Risk',      # OK
            'high risk',     # 제외
            'no call',       # 제외
            'Low Risk',      # OK
            'Low Risk',      # OK
            'Low Risk',      # OK
            'Low Risk',      # OK
            'Low Risk',      # OK
            'Low Risk',      # OK
            'Low Risk'       # OK
        ],
        'fetal_gender(gd_2)': ['XY'] * 15,
        'SeqFF': [10.0] * 15,
        'mapping_rate(%)': [98.0] * 15
    }
    
    df = pd.DataFrame(test_data)
    
    print("\n원본 데이터:")
    print(df[['sample_id', 'Result', 'MDResult']].to_string(index=False))
    print(f"\n총 샘플 수: {len(df)}")
    
    # 필터링 적용 (정확한 매칭)
    original_count = len(df)
    
    # Result 필터링 (High Risk) - exact match
    before = len(df)
    df = df[df['Result'].str.strip().str.lower() != 'high risk']
    print(f"\n1. Result 'High Risk' 제외: {before} -> {len(df)} ({before - len(df)} removed)")
    
    # Result 필터링 (No Call) - exact match
    before = len(df)
    df = df[df['Result'].str.strip().str.lower() != 'no call']
    print(f"2. Result 'No Call' 제외: {before} -> {len(df)} ({before - len(df)} removed)")
    
    # MDResult 필터링 (High Risk) - exact match
    before = len(df)
    mask = df['MDResult'].fillna('').str.strip().str.lower() != 'high risk'
    df = df[mask]
    print(f"3. MDResult 'High Risk' 제외: {before} -> {len(df)} ({before - len(df)} removed)")
    
    # MDResult 필터링 (No Call) - exact match
    before = len(df)
    mask = df['MDResult'].fillna('').str.strip().str.lower() != 'no call'
    df = df[mask]
    print(f"4. MDResult 'No Call' 제외: {before} -> {len(df)} ({before - len(df)} removed)")
    
    print("\n필터링 후 남은 데이터:")
    print(df[['sample_id', 'Result', 'MDResult']].to_string(index=False))
    print(f"\n최종 샘플 수: {len(df)} / {original_count}")
    
    # 검증
    # SAMPLE001, SAMPLE005, SAMPLE009, SAMPLE012 + SAMPLE013(None), SAMPLE014(Recalled), SAMPLE015(Normal)
    expected_count = 7
    
    print("\n" + "="*60)
    if len(df) == expected_count:
        print("✅ 테스트 성공!")
        print(f"   예상된 {expected_count}개 샘플만 남음")
        return True
    else:
        print("❌ 테스트 실패!")
        print(f"   예상: {expected_count}개, 실제: {len(df)}개")
        return False
    
if __name__ == "__main__":
    success = test_case_insensitive_filtering()
    sys.exit(0 if success else 1)

