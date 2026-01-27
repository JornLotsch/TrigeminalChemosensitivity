#!/bin/bash
# Batch convert all SVG files in current directory to 300 dpi PNG

for file in *.svg; do
    [ -e "$file" ] || continue   # skip if no svg files
    out="${file%.svg}.png"
    inkscape "$file" \
        --export-type=png \
        --export-filename="$out" \
        --export-dpi=300
    echo "Converted $file to $out (300 dpi)"
done

