# versions:
posseleff_simulations: commit 334008cd55c9d40c76fbea0318c4991e0e8d1fe1
ibdutils: commit 44e9838cd69afc94e1ccab8ef67869aa891e1f89


# running command 
nextflow ../main.nf -profile sge --num_reps 30 -resume


# files selected to download for plotting

fd -I '\.ne|member.pq' resdir/ > file_list.txt
fd -I '\.daf|.ibddist.ibdobj.gz|.true_ne' resdir/1000*rep0/ \
	>> file_list.txt
fd -I vcf.gz resdir/10003_sp_s03_rep0/vcf >> file_list.txt
fd -I '.ifm.ibdobj.gz|\.daf' \
	resdir/{20000_mp_s00_rep0,20001_mp_s01_rep0,20002_mp_s02_rep0,20003_mp_s03_rep0} \
	>> file_list.txt

## tar
tar -chzf posseleff_simu-r230420-sel.tgz -T file_list.txt

## rsync from local

rsync -a  lambda:/local/scratch/bing/posseleff_simulations/r230420/file_list.txt ./
rsync -azL lambda:/local/scratch/bing/posseleff_simulations/r230420 \
	--files-from file_list.txt ./
	