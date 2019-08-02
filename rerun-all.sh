#!/usr/bin/env bash

for i in out*; do for j in $i/seed*; do ./check.py $j --rerun; done; done
