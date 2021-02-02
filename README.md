# Hybrid NGMLR

Program requirement: 
minimap2, ngmlr, samtools and bamtools.

program usage:

    usage: python vulcan.py [-h] [-w WORK_DIR] [-t THREADS] [-i INPUT] [-o OUTPUT] [-p PERCENTILE] [-r REFERENCE] [-d] [-R] [-P]

    hybrid_NGMLR: a pipeline connects minimap2 and NGMLR

    optional arguments:
    -h, --help            show this help message and exit
    -w WORK_DIR, --WORK_DIR WORK_DIR
                            Directory of work directory
    -t THREADS, --threads THREADS
                            threads
    -i INPUT, --input INPUT
                            input read path
    -o OUTPUT, --output OUTPUT
                            output bamfile path
    -p PERCENTILE, --percentile PERCENTILE
                            percentile of cut-off, default: 50
    -r REFERENCE, --reference REFERENCE
                            raference path
    -d, --dry             only generate config
    -R, --raw_edit_distance
                            Use raw edit distance to do the cut-off
    -P, --pacbio          Input reads is pacbio reads

example: 

    python hybrid_ngmlr.py -w ../Saccharomyces_cerevisiae/work -r ../Saccharomyces_cerevisiae/GCF_000146045.2_R64_genomic.fna -i ../Saccharomyces_cerevisiae/reads_sv_simulated.fa -o ../Saccharomyces_cerevisiae/work/hybrid_ngmlr.bam -p 85 -t 60