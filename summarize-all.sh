#!/usr/bin/env bash

for i in out*; do for j in $i/seed*; do ./summarize.py $j; done; done
