#!/bin/bash

# Usage
if [ -z "$1" ]; then
    echo "Usage: $0 <max_samples_num>"
    exit 1
fi

max_samples=$1

# Detect total logical CPUs
total_cpus=$(lscpu | awk '/^CPU\(s\):/ {print $2}')

# Detect total system memory in GB
total_mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
total_mem_gb=$((total_mem_kb / 1024 / 1024))

# Reserve 20% memory for OS/other jobs
usable_mem_gb=$(echo "$total_mem_gb * 0.8" | bc)

# Allocate threads
bwa_threads=$(echo "$total_cpus * 0.85 / $max_samples" | bc)
samtools_threads=$(echo "$total_cpus * 0.15 / $max_samples" | bc)

# Enforce thread bounds
[ "$bwa_threads" -lt 4 ] && bwa_threads=4
[ "$samtools_threads" -lt 2 ] && samtools_threads=2
[ "$samtools_threads" -gt 8 ] && samtools_threads=8

# Calculate total RAM per sample
mem_per_sample=$(echo "$usable_mem_gb / $max_samples" | bc)

# Calculate memory per samtools thread, clamp to 1–4G
mem_per_thread=$(echo "$mem_per_sample / $samtools_threads" | bc)
[ "$mem_per_thread" -lt 1 ] && mem_per_thread=1
[ "$mem_per_thread" -gt 4 ] && mem_per_thread=4

# Output
echo "🧠  Detected $total_cpus logical CPUs"
echo "💾  Detected $total_mem_gb GB total memory"
echo "🔁  Targeting $max_samples samples in parallel"
echo ""
echo "✅  Recommended per-sample settings:"
echo "    BWA threads        : $bwa_threads"
echo "    Samtools threads   : $samtools_threads"
echo "    Samtools -m option : ${mem_per_thread}G"
echo ""
echo "⚠️  Total estimated samtools RAM: $((samtools_threads * mem_per_thread))G per sample × $max_samples samples = $((samtools_threads * mem_per_thread * max_samples))G"
