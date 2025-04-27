### multiBatch
f=GSE164378
echo "Run GSE164378"
python -u -W ignore ../scMAGCA/run_scMAGCA_batch.py --n_clusters 31 --n_batches 2 --data_file $f --alpha 1 --beta 1 --gama 0.01 --filter2
