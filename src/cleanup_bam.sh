#!/bin/bash

# cleanup_bam.sh - BAM 파일 정리 스크립트
# proper_paired.bam과 인덱스 파일만 남기고 나머지 BAM 파일들을 삭제

# 단일 샘플 디렉토리에서 BAM 파일 정리
cleanup_single_sample() {
    local sample_analysis_dir="$1"
    local sample_id="$2"
    local auto_confirm="$3"
    
    cd "$sample_analysis_dir" || return 1
    
    echo "Processing: $sample_analysis_dir"
    echo "Sample ID: $sample_id"
    
    # 삭제할 파일 패턴들 (proper_paired 제외)
    local files_to_delete=(
        "${sample_id}.sorted.bam"
        "${sample_id}.sorted.bam.bai"
        "${sample_id}.dedup.bam"
        "${sample_id}.dedup.bam.bai"
        "${sample_id}_dup.metrics"
        "${sample_id}.uniq.bam"
        "${sample_id}.uniq.bam.bai"
        "${sample_id}.of_orig.bam"
        "${sample_id}.of_orig.bam.bai"
        "${sample_id}.of_fetus.bam"
        "${sample_id}.of_fetus.bam.bai"
        "${sample_id}.of_mom.bam"
        "${sample_id}.of_mom.bam.bai"
        "${sample_id}.nf08_orig.bam"
        "${sample_id}.nf08_orig.bam.bai"
        "${sample_id}.nf08_fetus.bam"
        "${sample_id}.nf08_fetus.bam.bai"
        "${sample_id}.nf08_mom.bam"
        "${sample_id}.nf08_mom.bam.bai"
        "${sample_id}.nf09_orig.bam"
        "${sample_id}.nf09_orig.bam.bai"
        "${sample_id}.nf09_fetus.bam"
        "${sample_id}.nf09_fetus.bam.bai"
        "${sample_id}.nf09_mom.bam"
        "${sample_id}.nf09_mom.bam.bai"
    )
    
    local total_size=0
    local deleted_count=0
    local existing_files=()
    
    # 삭제 전 파일 크기 계산 및 존재하는 파일 목록 작성
    for file in "${files_to_delete[@]}"; do
        if [[ -f "$file" ]]; then
            existing_files+=("$file")
            # macOS와 Linux 호환성을 위한 파일 크기 확인
            if command -v stat >/dev/null 2>&1; then
                if [[ "$OSTYPE" == "darwin"* ]]; then
                    size=$(stat -f%z "$file" 2>/dev/null)
                else
                    size=$(stat -c%s "$file" 2>/dev/null)
                fi
                
                if [[ -n "$size" && "$size" != "0" ]]; then
                    total_size=$((total_size + size))
                    # numfmt가 있으면 사용, 없으면 바이트 단위로 표시
                    if command -v numfmt >/dev/null 2>&1; then
                        echo "  $file ($(numfmt --to=iec-i --suffix=B $size))"
                    else
                        echo "  $file ($size bytes)"
                    fi
                else
                    echo "  $file (size unknown)"
                fi
            else
                echo "  $file (size unknown)"
            fi
        fi
    done
    
    # 삭제할 파일이 없는 경우
    if [[ ${#existing_files[@]} -eq 0 ]]; then
        echo "  No BAM files to delete found."
        echo ""
        # 전역 변수로 결과 전달
        CLEANUP_DELETED_COUNT=0
        return 0
    fi
    
    # 파일 삭제
    for file in "${existing_files[@]}"; do
        if rm -f "$file" 2>/dev/null; then
            echo "  ✓ Deleted: $file"
            ((deleted_count++))
        else
            echo "  ✗ Failed to delete: $file"
        fi
    done
    
    echo "  Files deleted: $deleted_count"
    if command -v numfmt >/dev/null 2>&1 && [[ $total_size -gt 0 ]]; then
        echo "  Space freed: $(numfmt --to=iec-i --suffix=B $total_size)"
    else
        echo "  Space freed: $total_size bytes"
    fi
    
    # proper_paired 파일 확인
    if ls "${sample_id}.proper_paired.bam"* >/dev/null 2>&1; then
        echo "  ✓ Kept: $(ls ${sample_id}.proper_paired.bam*)"
    else
        echo "  ⚠ No proper_paired files found"
    fi
    
    echo ""
    # 전역 변수로 결과 전달
    CLEANUP_DELETED_COUNT=$deleted_count
    return 0
}

cleanup_bam_files() {
    local analysis_dir="$1"
    local sample_id="$2"
    local auto_confirm="${3:-false}"  # 자동 확인 모드 (기본값: false)
    
    # 전역 변수 초기화
    CLEANUP_DELETED_COUNT=0
    
    # 매개변수 검증
    if [[ -z "$analysis_dir" ]]; then
        echo "Usage: cleanup_bam_files <analysis_directory> [sample_id] [auto_confirm]"
        echo "Examples:"
        echo "  cleanup_bam_files /path/to/analysis 2505020003 true    # 특정 샘플"
        echo "  cleanup_bam_files /path/to/analysis                   # 모든 하위 디렉토리"
        return 1
    fi
    
    # 디렉토리 존재 확인
    if [[ ! -d "$analysis_dir" ]]; then
        echo "Error: Directory $analysis_dir does not exist"
        return 1
    fi
    
    echo "=========================================="
    echo "=== BAM Files Cleanup ==="
    echo "Analysis Directory: $analysis_dir"
    
    if [[ -n "$sample_id" ]]; then
        # 특정 샘플 처리
        local sample_analysis_dir="$analysis_dir/$sample_id"
        
        if [[ ! -d "$sample_analysis_dir" ]]; then
            echo "Error: Sample directory $sample_analysis_dir does not exist"
            return 1
        fi
        
        echo "Target: Single sample ($sample_id)"
        echo "Keeping: ${sample_id}.proper_paired.bam and its index file"
        echo "=========================================="
        
        # 단일 샘플에 대해서는 확인 프롬프트 표시
        if [[ "$auto_confirm" != "true" ]]; then
            echo "Files to be deleted:"
            # 미리 파일 목록 표시를 위해 임시로 single sample 함수 호출
            cd "$sample_analysis_dir" || return 1
            
            local files_to_delete=(
                "${sample_id}.sorted.bam" "${sample_id}.sorted.bam.bai"
                "${sample_id}.dedup.bam" "${sample_id}.dedup.bam.bai" "${sample_id}_dup.metrics"
                "${sample_id}.uniq.bam" "${sample_id}.uniq.bam.bai"
                "${sample_id}.of_orig.bam" "${sample_id}.of_orig.bam.bai"
                "${sample_id}.of_fetus.bam" "${sample_id}.of_fetus.bam.bai"
                "${sample_id}.of_mom.bam" "${sample_id}.of_mom.bam.bai"
                "${sample_id}.nf08_orig.bam" "${sample_id}.nf08_orig.bam.bai"
                "${sample_id}.nf08_fetus.bam" "${sample_id}.nf08_fetus.bam.bai"
                "${sample_id}.nf08_mom.bam" "${sample_id}.nf08_mom.bam.bai"
                "${sample_id}.nf09_orig.bam" "${sample_id}.nf09_orig.bam.bai"
                "${sample_id}.nf09_fetus.bam" "${sample_id}.nf09_fetus.bam.bai"
                "${sample_id}.nf09_mom.bam" "${sample_id}.nf09_mom.bam.bai"
            )
            
            local total_size=0
            for file in "${files_to_delete[@]}"; do
                if [[ -f "$file" ]]; then
                    if command -v stat >/dev/null 2>&1; then
                        if [[ "$OSTYPE" == "darwin"* ]]; then
                            size=$(stat -f%z "$file" 2>/dev/null)
                        else
                            size=$(stat -c%s "$file" 2>/dev/null)
                        fi
                        
                        if [[ -n "$size" && "$size" != "0" ]]; then
                            total_size=$((total_size + size))
                            if command -v numfmt >/dev/null 2>&1; then
                                echo "  $file ($(numfmt --to=iec-i --suffix=B $size))"
                            else
                                echo "  $file ($size bytes)"
                            fi
                        fi
                    else
                        echo "  $file"
                    fi
                fi
            done
            
            echo ""
            if command -v numfmt >/dev/null 2>&1 && [[ $total_size -gt 0 ]]; then
                echo "Total size to be freed: $(numfmt --to=iec-i --suffix=B $total_size)"
            else
                echo "Total size to be freed: $total_size bytes"
            fi
            echo ""
            
            read -p "Do you want to proceed with deletion? (y/N): " confirm
            if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
                echo "Operation cancelled."
                return 0
            fi
        else
            echo "Auto-confirm mode: Proceeding with deletion..."
        fi
        
        echo ""
        cleanup_single_sample "$sample_analysis_dir" "$sample_id" "$auto_confirm"
    else
        # 모든 하위 디렉토리 처리
        echo "Target: All subdirectories"
        echo "Scanning for sample directories..."
        echo "=========================================="
        
        local total_processed=0
        local total_deleted=0
        local processed_samples=()
        
        # analysis_dir의 하위 디렉토리들을 순회
        for subdir in "$analysis_dir"/*; do
            if [[ -d "$subdir" ]]; then
                local dir_name=$(basename "$subdir")
                
                # BAM 파일이 있는지 확인 (샘플 디렉토리인지 판단)
                if ls "$subdir"/*.bam >/dev/null 2>&1; then
                    echo "Found sample directory: $dir_name"
                    
                    # 디렉토리 이름을 샘플 ID로 사용
                    cleanup_single_sample "$subdir" "$dir_name" "true"
                    local cleanup_result=$?
                    
                    # cleanup_single_sample이 성공적으로 실행되었는지만 확인
                    # 삭제할 파일이 없어도 정상 처리된 것으로 간주
                    if [[ $cleanup_result -eq 0 ]]; then
                        processed_samples+=("$dir_name")
                        total_deleted=$((total_deleted + CLEANUP_DELETED_COUNT))
                        ((total_processed++))
                    else
                        echo "  ⚠ Failed to process directory: $dir_name"
                    fi
                else
                    echo "Skipping directory (no BAM files): $dir_name"
                fi
            fi
        done
        
        echo "=========================================="
        echo "=== Batch Cleanup Summary ==="
        echo "Total directories processed: $total_processed"
        echo "Total files deleted: $total_deleted"
        
        if [[ ${#processed_samples[@]} -gt 0 ]]; then
            echo "Processed samples:"
            for sample in "${processed_samples[@]}"; do
                echo "  - $sample"
            done
        else
            echo "No sample directories found with deletable BAM files"
        fi
    fi
    
    echo "=========================================="
    echo "Cleanup completed!"
    return 0
}

# 스크립트가 직접 실행된 경우
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # 명령행 인자 처리
    if [[ $# -lt 1 ]]; then
        echo "Usage: $0 <analysis_directory> [sample_id] [auto_confirm]"
        echo "Examples:"
        echo "  $0 /home/ken/ken-nipt/analysis/250430_01 2505020003 true    # 특정 샘플"
        echo "  $0 /home/ken/ken-nipt/analysis/250430_01                   # 해당 디렉토리의 모든 하위 샘플"
        echo "  $0 /home/ken/ken-nipt/analysis                             # analysis 디렉토리의 모든 work 디렉토리"
        exit 1
    fi
    
    cleanup_bam_files "$1" "$2" "$3"
fi
