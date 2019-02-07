#!/bin/bash

col1=$1
col2=$2
col3=$3

sort -n -k $col1 | sort -s -n -k $col2 | sort -s -n -k $col3 $4
