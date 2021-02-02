# Vulcan

Program requirement: 
minimap2, ngmlr, samtools=1.9 and bamtools.

program usage:

    usage: vulcan.py [-h] -i INPUT [INPUT ...] -r REFERENCE -o OUTPUT [-w WORK_DIR] [-t THREADS] [-p PERCENTILE] [-f] [-d] [-R]
                    (-clr | -hifi | -ont | -any | -hclr | -hhifi | -hont | -cmd)

    Vulcan: Map long and prosper🖖, a pipeline melds minimap2 and NGMLR

    optional arguments:
    -h, --help            show this help message and exit
    -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                            input read path, can accept multiple files
    -r REFERENCE, --reference REFERENCE
                            reference path
    -o OUTPUT, --output OUTPUT
                            vulcan's output bamfile path
    -w WORK_DIR, --work_dir WORK_DIR
                            Directory of work, store temp files, default: ./vulcan_work
    -t THREADS, --threads THREADS
                            threads, default: 1
    -p PERCENTILE, --percentile PERCENTILE
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

example: 

    python vulcan.py -r {reference} -i {input} -o {output} -t 60 -any