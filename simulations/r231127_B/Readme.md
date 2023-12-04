# Goal

simulation inbreeding by changing selfing rate (single population model)


# Command

```
python3 generate_json.py
nextflow ../../main.nf \
    -profile sge  \
    -entry WF_SP \
    -resume \
    --sp_sets_json sp_genome_sets.json \
    --ibdne_no_diploid_convertion true \
    --peak_validate_meth ihs \
    --num_reps 1 
```

