#!/bin/bash

PAF=$1

sort -V -k1,1 -k6,6 $PAF | sponge $PAF
