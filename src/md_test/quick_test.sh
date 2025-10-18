#!/usr/bin/env bash
# Quick test: Generate 2 samples to verify everything works
set -euo pipefail

cd "$(dirname "$0")"

echo "========================================="
echo "Quick Test: 2 Artificial Samples"
echo "========================================="
echo ""
echo "Settings:"
echo "  Disease: 1p36 deletion syndrome"
echo "  FF: 5%"
echo "  Coverage: 10M pairs"
echo "  Samples: 2 moms × 1 fetus = 2 samples"
echo "  Expected time: ~10 minutes"
echo ""
echo "Press Ctrl+C to cancel, or wait 3 seconds to start..."
sleep 3
echo ""

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Clean previous test output
rm -rf "$SCRIPT_DIR/test_output" "$SCRIPT_DIR/test_logs" "$SCRIPT_DIR/summary"

# Run minimal test using TSV files (new interface)
python3 "$SCRIPT_DIR/run_parallel_v2.py" \
  --mom_tsv "$SCRIPT_DIR/mom_list.tsv" \
  --female_tsv "$SCRIPT_DIR/female_fetus_list.tsv" \
  --male_tsv "$SCRIPT_DIR/male_fetus_list.tsv" \
  --md_bed "$SCRIPT_DIR/test_1p36_only.bed" \
  --script "$SCRIPT_DIR/make_artificial.sh" \
  --output "$SCRIPT_DIR/test_output" \
  --ff_targets 5 \
  --coverages 10M \
  --n_moms 2 \
  --n_fetuses 1 \
  --workers 1

echo ""
echo "========================================="
echo "Test Complete!"
echo "========================================="
echo ""

# Show results
if [ -d "$SCRIPT_DIR/test_output" ]; then
  echo "✓ Generated files:"
  find "$SCRIPT_DIR/test_output" -name "*.bam" | sort | while read f; do
    size=$(du -h "$f" | cut -f1)
    echo "  $f ($size)"
  done
  echo ""
  
  count=$(find "$SCRIPT_DIR/test_output" -name "*.bam" | wc -l)
  echo "✓ Total BAM files: $count"
  echo ""
  
  if [ -f "$SCRIPT_DIR/summary/all_samples.tsv" ]; then
    echo "✓ Sample mapping:"
    column -t -s $'\t' "$SCRIPT_DIR/summary/all_samples.tsv"
    echo ""
  fi
  
  if [ -f "$SCRIPT_DIR/test_logs/1p36/FF05_10M/"*.log ]; then
    echo "✓ Log files:"
    ls -1 "$SCRIPT_DIR/test_logs/1p36/FF05_10M/"*.log
    echo ""
  fi
  
  echo "Next steps:"
  echo "  1. Verify BAM files: samtools view $SCRIPT_DIR/test_output/*/FF*/*.bam | head"
  echo "  2. Check logs: less $SCRIPT_DIR/test_logs/*/FF*/*.log"
  echo "  3. Run larger test with more samples"
else
  echo "✗ No output directory found. Check errors above."
  exit 1
fi

