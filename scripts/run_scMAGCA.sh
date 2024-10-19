### CITE-Seq 
f=10x1kpbmc
echo "Run 10x1kpbmc"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters 5 --data_file $f --save_dir 10x1kpbmc --filter2
f=10x5kpbmc
echo "Run 10x5kpbmc"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters 12 --data_file $f --save_dir 10x5kpbmc --filter2
f=10x5kpbmcTotalSeq
echo "Run 10x5kpbmcTotalSeq"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters 12 --data_file $f --save_dir 10x5kpbmcTotalSeq --filter2
f=10xmalt
echo "Run 10xmalt"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters 11 --data_file $f --save_dir 10xmalt --filter2
f=GSE128639_BMNC
echo "Run GSE128639_BMNC"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters 27 --data_file $f --save_dir GSE128639_BMNC --filter2
f=spector
echo "Run spector"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters 16 --data_file $f --save_dir spector --filter2
### multiBatch
f=GSE164378
echo "Run GSE164378"
python -u -W ignore ../src/run_scMAGCA_batch.py --n_clusters 32 --data_file $f --save_dir GSE164378 --filter2



### SMAGE-Seq
f=human_brain_3k
echo "Run human_brain_3k"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters 10 --data_file $f --save_dir human_brain_3k --filter1 --filter2 --ad_out 32
f=human_pbmc_3k
echo "Run human_pbmc_3k"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters 8 --data_file $f --save_dir human_pbmc_3k --filter1 --filter2 --ad_out 32
f=mouse_brain_5k
echo "Run mouse_brain_5k"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters 11 --data_file $f --save_dir mouse_brain_5k --filter1 --filter2 --ad_out 32
### no-label
f=GSM4949911_tea
echo "Run GSM4949911_tea"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters -1 --data_file $f --save_dir GSM4949911_tea --filter1 --filter2 --ad_out 32
f=pbmc_10x_public
echo "pbmc_10x_public"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters -1 --data_file $f --save_dir pbmc_10x_public --filter1 --filter2 --ad_out 32
f=10x-Multiome-Pbmc10k
echo "Run 10x-Multiome-Pbmc10k"
python -u -W ignore ../src/run_scMAGCA.py --n_clusters -1 --data_file $f --save_dir 10x-Multiome-Pbmc10k --filter1 --filter2 --ad_out 32


#f=Kindy
#echo "Run Kindy"
#python -u -W ignore ../src/run_scGMD.py --n_clusters -1 --ae_weight_file AE_weights_Muto-2021.pth.tar --data_file $f --save_dir Kindy --embedding_file --prediction_file --filter1 --filter2 --ad_out 32
#f=MouseKindy
#echo "Run MouseKindy"
#python -u -W ignore ../src/run_scGMD.py --n_clusters -1 --ae_weight_file AE_weights_Muto-2021.pth.tar --data_file $f --save_dir MouseKindy --embedding_file --prediction_file --filter1 --filter2 --f1 2000 --f2 2000

