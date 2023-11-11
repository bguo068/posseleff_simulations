# Goal

- incoporate inbreeding to multiple simulation model
- allow to simulate only chromosome for a single chromosome

# command

```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -profile sge \
    -entry WF_MP \
    -resume \
    --ibdne_no_diploid_convertion true \
    --num_reps 1 \
    --sim_inbreeding true 
```

# finding

Single chrome genome simulation seems to capture more fine-scale structure
instead of the larger true population structure. This is likely due to a lack
of average-out effects from multiple-chromsome genome simulations or multiple
independent chromosome simulations.
