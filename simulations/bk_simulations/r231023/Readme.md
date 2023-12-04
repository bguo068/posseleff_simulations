# main changes compared to `inbreed` branch

>
> - conditional on the non fixing mutation to have allele frequency > 0.2
> - change default `max_restart ` from 100 to 100000 within the `./bin/sim_single_pop_whole_genome.py`

in addition to r231020, the time span of selection is reduced from the default 80 to 50, as
80 might be to long so that the selection signal is decreasing

# Command
```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -entry WF_SP \
    -resume \
    -profile sge \
    --ibdne_no_diploid_convertion true \
    --sim_inbreeding true
```
