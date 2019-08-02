#!/usr/bin/env bash

for i in out-sp2-gen4-loci100*; do for j in $i/seed*; do ./summarize.py $j; done; done
