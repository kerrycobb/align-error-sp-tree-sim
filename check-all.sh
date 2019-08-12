#!/usr/bin/env bash

for i in out*; do for j in $i/seed*; do echo $j; ./check.py $j; done; done
