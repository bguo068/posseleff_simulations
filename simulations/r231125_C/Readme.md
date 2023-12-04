# Goal

inbreeding via selfing (Multiple population model)


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
