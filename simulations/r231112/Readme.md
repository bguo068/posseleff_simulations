# Goal

- simulate each chromosome independently
- explore assortative mating via `modifyChild` instead of `mateChoice`

# edits

main.nf
single_pop.slim



# command

```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -entry WF_SP \
    -resume \
    -profile sge \
    --ibdne_no_diploid_convertion true \
    --num_reps 1 

```

# findings and notes

1. selection correction is indeed inbreeding-dependent
    - selection bias in simulation with lower amount of inbreeding is correctable
    - selection bias in simulation with higher amount of inbreeding is harder to correct and can cause overcorrection

2. repeat might be needed to see the general trend of selection effect and bias correction effects
    - especially for the non-inbreeding simulation cases
     