### multiBatch
f=GSE164378
echo "Run GSE164378"
python -u -W ignore ../src/run_scMAGCA_batch.py --n_clusters 32 --data_file $f --save_dir GSE164378 --filter2