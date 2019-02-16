#!/bin/bash

xterm -e "diff $1 $2 |& less"
