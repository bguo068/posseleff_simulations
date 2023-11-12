# Goal

relative to r231110:
- return to simulating 14 chromosomes independently under multiple-population model
- use modifyChild instead of mateChoice as the former is faster
    - Also, modifyChild can avoid extreme assortative mating-based inbreeding 
- search over multiple parameters

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

# genome sets

- json: `./genome_sets.json`
- groovy script to generate the json file: `mpsets.groovy`

# finding

Got expected results:
- non-inbreeding: removing peaks improves pop structure inference
- weak-inbreeding: removing peaks improves pop structure inference
- strong-inbreeding: removing peaks does not improve pop structure inference
