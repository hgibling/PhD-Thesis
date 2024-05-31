set -e
true
true
/.mounts/labs/awadallalab/private/hgibling/SOFT/SPAdes/3.14.1/SPAdes-3.14.1-Linux/bin/spades-hammer /.mounts/labs/simpsonlab/users/hgibling/graphs/spades-test/corrected/configs/config.info
/.mounts/labs/simpsonlab/sw/miniconda3/bin/python /.mounts/labs/awadallalab/private/hgibling/SOFT/SPAdes/3.14.1/SPAdes-3.14.1-Linux/share/spades/spades_pipeline/scripts/compress_all.py --input_file /.mounts/labs/simpsonlab/users/hgibling/graphs/spades-test/corrected/corrected.yaml --ext_python_modules_home /.mounts/labs/awadallalab/private/hgibling/SOFT/SPAdes/3.14.1/SPAdes-3.14.1-Linux/share/spades --max_threads 16 --output_dir /.mounts/labs/simpsonlab/users/hgibling/graphs/spades-test/corrected --gzip_output
true
/.mounts/labs/simpsonlab/sw/miniconda3/bin/python /.mounts/labs/awadallalab/private/hgibling/SOFT/SPAdes/3.14.1/SPAdes-3.14.1-Linux/share/spades/spades_pipeline/scripts/breaking_scaffolds_script.py --result_scaffolds_filename /.mounts/labs/simpsonlab/users/hgibling/graphs/spades-test/scaffolds.fasta --misc_dir /.mounts/labs/simpsonlab/users/hgibling/graphs/spades-test/misc --threshold_for_breaking_scaffolds 3
true
