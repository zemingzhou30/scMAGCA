### multiBatch
f=GSE164378
echo "Run GSE164378"
python -u -W ignore ../src/run_scMAGCA_batch.py --n_clusters 31 --data_file $f --filter2
