### CITE-Seq 
f=10x1kpbmc
echo "Run 10x1kpbmc"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --n_clusters 5 --data_file $f --alpha 0.5 --beta 1 --gama 0.1 --filter2
f=10x5kpbmc
echo "Run 10x5kpbmc"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --n_clusters 12 --data_file $f --alpha 1 --beta 0.5 --gama 0.01 --filter2
f=10x5kpbmcTotalSeq
echo "Run 10x5kpbmcTotalSeq"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --n_clusters 12 --data_file $f --alpha 1 --beta 0.5 --gama 0.01 --filter2
f=10xmalt
echo "Run 10xmalt"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --n_clusters 11 --data_file $f --alpha 0.5 --beta 1 --gama 0.05 --filter2
f=GSE128639_BMNC
echo "Run GSE128639_BMNC"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --n_clusters 27 --data_file $f --alpha 1 --beta 0.5 --gama 0.01 --filter2
f=spector
echo "Run spector"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --n_clusters 16 --data_file $f --alpha 0.5 --beta 1 --gama 0.02 --filter2



### SMAGE-Seq
f=human_brain_3k
echo "Run human_brain_3k"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --n_clusters 10 --data_file $f --alpha 1 --beta 1 --gama 0.2 --filter1 --filter2 --ad_out 32
f=human_pbmc_3k
echo "Run human_pbmc_3k"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --n_clusters 8 --data_file $f --alpha 1 --beta 1 --gama 0.3 --filter1 --filter2 --ad_out 32
f=mouse_brain_5k
echo "Run mouse_brain_5k"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --n_clusters 11 --data_file $f --alpha 1 --beta 1 --gama 0.2 --filter1 --filter2 --ad_out 32
### no-label
f=GSM4949911_tea
echo "Run GSM4949911_tea"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --data_file $f --filter1 --filter2 --ad_out 32
f=pbmc_10x_public
echo "pbmc_10x_public"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --data_file $f --filter1 --filter2 --ad_out 32
f=10x-Multiome-Pbmc10k
echo "Run 10x-Multiome-Pbmc10k"
python -u -W ignore ../scMAGCA/run_scMAGCA.py --data_file $f --filter1 --filter2 --ad_out 32
