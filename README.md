# Vulcan

[![Anaconda-Server Badge](https://anaconda.org/bioconda/vulcan/badges/version.svg)](https://anaconda.org/bioconda/vulcan)  [![Anaconda-Server Badge](https://anaconda.org/bioconda/vulcan/badges/downloads.svg)](https://anaconda.org/bioconda/vulcan)

Installation: 

    conda install -c bioconda vulcan

program usage:
    
    usage: vulcan [-h] -i INPUT [INPUT ...] -r REFERENCE -o OUTPUT [-w WORK_DIR] [-t THREADS] [-p PERCENTILE [PERCENTILE ...]] [-f] [-d] [-R] [-clr | -hifi | -ont | -any | -hclr | -hhifi | -hont | -cmd]


    vulcan: map long reads and prosper🖖, a long read mapping pipeline that melds minimap2 and NGMLR
    
    optional arguments:
    -h, --help            show this help message and exit
    -w WORK_DIR, --work_dir WORK_DIR
                            Directory of work, store temp files, default: ./vulcan_work
    -t THREADS, --threads THREADS
                            threads, default: 1
    -p PERCENTILE [PERCENTILE ...], --percentile PERCENTILE [PERCENTILE ...]
                            percentile of cut-off, default: 90
    -f, --full            keep all temp file
    -d, --dry             only generate config
    -R, --raw_edit_distance
                            Use raw edit distance to do the cut-off
    -clr, --pacbio_clr    Input reads is pacbio CLR reads
    -hifi, --pacbio_hifi  Input reads is pacbio hifi reads
    -ont, --nanopore      Input reads is Nanopore reads
    -any, --anylongread   Don't know which kind of long read
    -hclr, --humanclr     Human pacbio CLR read
    -hhifi, --humanhifi   Human pacbio hifi reads
    -hont, --humannanopore
                            Human Nanopore reads
    -cmd, --custom_cmd    Use minimap2 and NGMLR with user's own parameter setting

    Required arguments::
    -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                            input read path, can accept multiple files
    -r REFERENCE, --reference REFERENCE
                            reference path
    -o OUTPUT, --output OUTPUT
                            vulcan's output's prefix, the output will be prefix_{percentile}.bam

example: 

    ./vulcan -r ./test/GCF_000146045.2_R64_genomic.fna -i ./test/test_reads.fa -w ./test/ -o vulcan
