#!/bin/bash

set -e

for tex_path in ./grid*.tex
do
    pdf_path="${tex_path%.*}.pdf"
    cropped_pdf_path="$(dirname $pdf_path)/cropped-$(basename $pdf_path)"
    latexmk -C "$tex_path"
    latexmk -pdf "$tex_path"
    pdfcrop "$pdf_path" "$cropped_pdf_path"
    latexmk -C "$tex_path"
done
