# Notes on Simulations

## Subfolders

1. **Simulations Without Explicit Inbreeding Modeling**
   - `r231116`

2. **Simulations With Explicit Inbreeding Modeling**
   - **Modeling Inbreeding via Shrinking Population Size**
     - Multiple-Population Model: `r231125_A`
     - Single-Population Model: `r231127_A`
   - **Modeling Inbreeding via Positive Assortative Mating**
     - Multiple-Population Model: `r231125_B2`
     - Single-Population Model: `r231127_C`
   - **Modeling Inbreeding via Selfing**
     - Multiple-Population Model: `r231125_C`
     - Single-Population Model: `r231127_B`

## Parameters Used for Each Set of Simulations

Simulation parameters are configured at three levels:

1. **Global Defaults**: 
   - Defined in the `main.nf` file as `sp_defaults` and `mp_defaults`. 
   - These set default values for all genomes in single/multiple population models and are used as a base in the subsequent levels.

2. **Within-Pipeline Parameter Sets**:
   - Named as `sp_sets` and `mp_sets`.
   - These are "hardcoded" simulation parameters that build upon the global defaults.

3. **Custom JSON Files**:
   - When JSON files are provided, they override the `sp_sets` and `mp_sets`.
   - If parameters not specified in these files, nextflow will access `sp_defaults` and `mp_defaults` for default values.

### Specific Notes:

- For the `r231116` subfolder, parameters are specified at the second level. Refer to `main.nf` for detailed information.
- For the remaining subfolders (`r231125_*` and `r231127_*`), simulation parameters are provided in JSON files located within each subfolder, such as `r231125_A/mp_genome_sets.json` and `r231127_A/sp_genome_sets.json`.

### Parameter names
Some parameters used in the simulation script, espeically for the assortative mating model, have different names in the Supplementary Note 2. Here is a map of them.

| In Paper        | in script                        |
| --------------- | -------------------------------- |
| G               |sim_relatedness_g                 |
| B               |sim_relatedness_bypass            |
| C               |sim_relatedness_bypass_complement |
| D               |sim_relatedness_power             |
| \delta          |sim_relatedness_delta             |
