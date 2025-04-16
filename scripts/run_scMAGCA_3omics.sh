### three-omics 
f=GSE158013
echo "Run GSE158013"
python -u -W ignore ../src/run_scMAGCA_3omics.py --n_clusters -1 --data_file $f --save_dir GSE158013 --embedding_file --prediction_file --filter2 --filter3 --f2 2000 --f3 2000 --ad_out 32