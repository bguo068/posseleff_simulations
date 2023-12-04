# main changes compared to main (554fce4b6a4ff6b87b05a9c6e75ad71251fb5da0)
- using `sim_single_pop_whole_genome.py`
- run single population only for this time
- focusing how inbreading starting time, power and delta affect the results


# Command
```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -entry WF_SP \
    -resume \
    -profile sge \
    --ibdne_no_diploid_convertion true \
    --sim_inbreeding true
```
# Process of a few simulation sets failed with no sucess
```
WF_SP:PROC_DIST_NE (140002)
WF_SP:PROC_DIST_NE (140004)
WF_SP:PROC_DIST_NE (240000)
WF_SP:RUN_IBDNE (140000_false)
WF_SP:RUN_IBDNE (340000_false)
```

some files are not copied due to interupted process are manually copied
```
cat /local/scratch/bing/posseleff_simulations/simulations/r231014/resdir/pipeline_trace.txt | grep PROC_DIST_NE | grep 140002 | head -1
rsync -aL  \
    /autofs/scratch/bing/posseleff_simulations/simulations/r231014/work/b6/a3565b1f55892203de34c3a252757b/140002_2.0_10.0_none.ibddist.ibdobj.gz \
    /local/chib/toconnor_grp/bing/posseleff_simulations/simulations/r231014/resdir/140002_sp_rels03_14b/ibddist_ibd/
cat /local/scratch/bing/posseleff_simulations/simulations/r231014/resdir/pipeline_trace.txt | grep PROC_DIST_NE | grep 140004 | head -1
rsync -aL /autofs/scratch/bing/posseleff_simulations/simulations/r231014/work/95/71df82f7210f63a1fcd8be21223913/140004_2.0_10.0_none.ibddist.ibdobj.gz /local/chib/toconnor_grp/bing/posseleff_simulations/simulations/r231014/resdir/140004_sp_rels03_14c/ibddist_ibd/
cat /local/scratch/bing/posseleff_simulations/simulations/r231014/resdir/pipeline_trace.txt | grep PROC_DIST_NE | grep 240000 | head -1
rsync -aL /autofs/scratch/bing/posseleff_simulations/simulations/r231014/work/8f/77836e4528a3a785e45dc63b898621/240000_2.0_10.0_none.ibddist.ibdobj.gz  /local/chib/toconnor_grp/bing/posseleff_simulations/simulations/r231014/resdir/240000_sp_rels03_24a/ibddist_ibd/
```