#!/bin/bash

search_dir="workflow/reports/paper_figures/"
output_dir="workflow/reports/paper_figures_pdf/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

#copy over the pdfs
for infile in "$search_dir"*.pdf
do
  outfile="$output_dir"$(basename "$infile")
  echo "Copying ${infile} to ${outfile}"
  cp "$infile" "$outfile"
done

#convert the svgs
for infile in "$search_dir"*.svg
do
  outfile="$output_dir"$(basename "$infile" .svg).pdf
  echo "Converting ${infile} to ${outfile}"
  cat $infile | inkscape --pipe --export-filename=$outfile -d 600
done