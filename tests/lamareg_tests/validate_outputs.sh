#!/bin/bash
#
# Validate LAMAReg Registration Outputs
# Checks existing registration outputs for completeness and quality
#
# Usage: ./validate_outputs.sh <registration_output_dir>
#

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}LAMAReg Output Validation${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

if [ $# -lt 1 ]; then
    echo -e "${RED}Error: Registration output directory required${NC}"
    echo "Usage: $0 <registration_output_dir>"
    echo ""
    echo "This script checks for:"
    echo "  - Presence of all required output files"
    echo "  - File integrity (non-zero size, valid headers)"
    echo "  - DICE scores from QC CSV files"
    echo "  - Transformation file validity"
    exit 1
fi

OUTPUT_DIR="$1"
VALIDATION_REPORT="$OUTPUT_DIR/validation_report.txt"

# Initialize report
echo "LAMAReg Output Validation Report" > "$VALIDATION_REPORT"
echo "Validation Date: $(date)" >> "$VALIDATION_REPORT"
echo "Output Directory: $OUTPUT_DIR" >> "$VALIDATION_REPORT"
echo "========================================" >> "$VALIDATION_REPORT"
echo "" >> "$VALIDATION_REPORT"

# Counters
FILES_CHECKED=0
FILES_VALID=0
FILES_INVALID=0
FILES_MISSING=0

# Required file patterns
REQUIRED_PATTERNS=(
    "*0GenericAffine.mat"
    "*1Warp.nii.gz"
    "*1InverseWarp.nii.gz"
    "*2Warp.nii.gz"
    "*2InverseWarp.nii.gz"
    "*_fixed_parc.nii.gz"
    "*_moving_parc.nii.gz"
    "*_registered_parc.nii.gz"
    "*_dice_scores.csv"
    "*Warped.nii.gz"
)

# Check function
check_file_pattern() {
    local pattern="$1"
    local files=("$OUTPUT_DIR"/$pattern)
    
    ((FILES_CHECKED++))
    
    if [ ! -e "${files[0]}" ]; then
        echo -e "${YELLOW}⚠ MISSING${NC}: $pattern"
        echo "⚠ MISSING: $pattern" >> "$VALIDATION_REPORT"
        ((FILES_MISSING++))
        return 1
    fi
    
    for file in "${files[@]}"; do
        if [ -f "$file" ]; then
            local filesize=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null)
            local filename=$(basename "$file")
            
            if [ "$filesize" -lt 100 ]; then
                echo -e "${RED}✗ INVALID${NC}: $filename (size: $filesize bytes)"
                echo "✗ INVALID: $filename (size: $filesize bytes)" >> "$VALIDATION_REPORT"
                ((FILES_INVALID++))
            else
                echo -e "${GREEN}✓ VALID${NC}: $filename ($(numfmt --to=iec-i --suffix=B $filesize 2>/dev/null || echo "${filesize} bytes"))"
                echo "✓ VALID: $filename" >> "$VALIDATION_REPORT"
                ((FILES_VALID++))
                
                # Additional validation for NIfTI files
                if [[ "$file" == *.nii.gz ]] && command -v fslinfo &>/dev/null; then
                    if ! fslinfo "$file" &>/dev/null; then
                        echo -e "  ${YELLOW}  Warning: Invalid NIfTI header${NC}"
                        echo "    Warning: Invalid NIfTI header" >> "$VALIDATION_REPORT"
                    fi
                fi
            fi
        fi
    done
}

# Check DICE scores
check_dice_scores() {
    local csv_files=("$OUTPUT_DIR"/*_dice_scores.csv)
    
    if [ ! -e "${csv_files[0]}" ]; then
        echo -e "${YELLOW}⚠ No DICE score CSV files found${NC}"
        return
    fi
    
    echo ""
    echo -e "${BLUE}DICE Score Summary:${NC}"
    echo "" >> "$VALIDATION_REPORT"
    echo "DICE Score Summary:" >> "$VALIDATION_REPORT"
    
    for csv_file in "${csv_files[@]}"; do
        if [ -f "$csv_file" ]; then
            local filename=$(basename "$csv_file")
            echo -e "${BLUE}  $filename:${NC}"
            echo "  $filename:" >> "$VALIDATION_REPORT"
            
            if command -v python3 &>/dev/null; then
                python3 << EOF | tee -a "$VALIDATION_REPORT"
import csv
import sys

try:
    with open('$csv_file', 'r') as f:
        reader = csv.DictReader(f)
        scores = []
        for row in reader:
            if 'dice' in row:
                dice = float(row['dice'])
                label = row.get('label', row.get('region', 'unknown'))
                scores.append(dice)
                print(f"    {label}: {dice:.3f}")
        
        if scores:
            avg = sum(scores) / len(scores)
            print(f"    Average: {avg:.3f}")
            if avg >= 0.70:
                print(f"    Status: GOOD (>= 0.70)")
            elif avg >= 0.60:
                print(f"    Status: ACCEPTABLE (>= 0.60)")
            else:
                print(f"    Status: POOR (< 0.60)")
except Exception as e:
    print(f"    Error reading CSV: {e}")
EOF
            else
                echo "    (Python3 not available for detailed analysis)"
                echo "    Line count: $(wc -l < "$csv_file")"
            fi
            echo ""
        fi
    done
}

# Run validation
echo "Checking for required output files..."
echo ""

for pattern in "${REQUIRED_PATTERNS[@]}"; do
    check_file_pattern "$pattern"
done

# Check DICE scores
check_dice_scores

# Summary
echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Validation Summary${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "Files checked: ${FILES_CHECKED}"
echo -e "${GREEN}Valid: ${FILES_VALID}${NC}"
echo -e "${RED}Invalid: ${FILES_INVALID}${NC}"
echo -e "${YELLOW}Missing: ${FILES_MISSING}${NC}"
echo ""

echo "" >> "$VALIDATION_REPORT"
echo "========================================" >> "$VALIDATION_REPORT"
echo "Validation Summary:" >> "$VALIDATION_REPORT"
echo "  Files checked: $FILES_CHECKED" >> "$VALIDATION_REPORT"
echo "  Valid: $FILES_VALID" >> "$VALIDATION_REPORT"
echo "  Invalid: $FILES_INVALID" >> "$VALIDATION_REPORT"
echo "  Missing: $FILES_MISSING" >> "$VALIDATION_REPORT"
echo "" >> "$VALIDATION_REPORT"

echo "Detailed report: $VALIDATION_REPORT"
echo ""

if [ $FILES_INVALID -eq 0 ] && [ $FILES_MISSING -eq 0 ]; then
    echo -e "${GREEN}✓ All outputs are valid!${NC}"
    echo "Result: ALL OUTPUTS VALID" >> "$VALIDATION_REPORT"
    exit 0
else
    echo -e "${RED}✗ Some outputs are invalid or missing${NC}"
    echo "Result: VALIDATION FAILED" >> "$VALIDATION_REPORT"
    exit 1
fi
