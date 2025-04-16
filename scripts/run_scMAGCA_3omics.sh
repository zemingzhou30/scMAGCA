### three-omics 
f=GSE158013
echo "Run GSE158013"
python -u -W ignore ../src/run_scMAGCA_3omics.py --data_file $f --filter2 --filter3 --f2 2000 --f3 2000 --ad_out 32