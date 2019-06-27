#!/bin/bash

colNum=$1

awk -v "col=$1" '{sum += $col} END {OFMT="%20.15e";print sum}' 
