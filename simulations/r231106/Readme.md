# main changes compared to `inbreed2` branch

>
> - conditional on the non fixing mutation to have allele frequency > 0.2
> - change default `max_restart ` from 100 to 100000 within the `./bin/sim_single_pop_whole_genome.py`

> in addition to r231020, the time span of selection is reduced from the default 80 to 50, as
> 80 might be to long so that the selection signal is decreasing

in addition to r231023
- change the delta and inbreeding time g
- allow excluding chromosome with unestablished for the Ne estimation purpose

# commits
- slim: 24c544fdebe90333541fcf48a988c3ebe8c7f2f2
- ibdutils: 0c4a339ddbe19fe804dd54cb3e41e3e4685def0d
- posseleff_simulations: commit 85b9a71a0fea8395b0d7d0e5cb056ad5b6b10dd6 (HEAD -> inbreeding2)


# Command
```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -entry WF_SP \
    -resume \
    -profile sge \
    --ibdne_no_diploid_convertion true \
    --sim_inbreeding true
```
