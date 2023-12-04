# Goal

simulate high inbreeding by decreasing population size for multiple-pop model

# Command

```
python3 generate_json.py
nextflow ../../main.nf \
    -profile sge \
    -entry WF_MP \
    -resume \
    --mp_sets_json mp_genome_sets.json \
    --ibdne_no_diploid_convertion true \
    --peak_validate_meth ihs \
    --num_reps 1 \
    --resdir resdir
```

# Note:

The subset of simulations with sel_mig = 0.01 were in downstream analyses
