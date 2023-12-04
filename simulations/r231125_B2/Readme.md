# Goal

simulate high inbreeding by assortative mating for multiple pop model

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

Non-default bypass parameters were used to allow higher inbreeding