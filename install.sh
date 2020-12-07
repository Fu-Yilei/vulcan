#!/bin/bash
conda create --name hybrid_NGMLR
conda activate hybrid_NGMLR

conda install ngmlr minimap2 SURVIVOR
conda install  samtools=1.9

echo "Install complete"