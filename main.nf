#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

params.ifm_transform = ["square", "cube", "none"][0]
params.ifm_ntrials = 1000
params.test = false
params.num_reps = 1
params.resdir = "res"

params.ibdne_mincm = 2
params.ibdne_minregion = 10
params.ibdne_flatmeth = 'none'
// default to false so it is consistent with previous version
// can be set to true on the nextflow command line
params.ibdne_no_diploid_convertion = false
params.sim_inbreeding = false

def resdir = params.resdir

def sp_defaults = [
    seqlen : 100 * 15000,
    selpos : Math.round(0.33 * 100 * 15000),
    num_origins : 1,
    N : 10000,
    h : 0.5,
    s : 0.3,
    g_sel_start : 80,
    r : 0.01 / 15_000,
    sim_relatedness : 0,
    g_ne_change_start : 200,
    N0 : 1000,
    u : 1e-8,
    nsam : 1000, // haploid
]

def mp_defaults = [
    seqlen : 100 * 15000,
    selpos : Math.round(100 * 15000 * 0.33),
    num_origins : 1,
    N : 10000,
    h : 0.5,
    s : 0.3,
    g_sel_start : 80,
    r : 0.01 / 15_000,
    sim_relatedness : 0,
    mig : 1e-5,
    sel_mig : 0.01,
    npop : 5,
    nsam : 200, // haploid
    Tsplit : 500,
    u : 1e-8,
]

def sp_sets = [
    sp_neu: sp_defaults + [s: 0.0, genome_set_id: 10000],
    // sp_s01: sp_defaults + [s: 0.1, genome_set_id: 10001],
    // sp_s02: sp_defaults + [s: 0.2, genome_set_id: 10002],
    sp_s03: sp_defaults + [s: 0.3, genome_set_id: 10003],
    // sp_g040:sp_defaults + [g_sel_start: 40, genome_set_id: 10004],
    // sp_g080:sp_defaults + [g_sel_start: 80, genome_set_id: 10005],
    // sp_g120:sp_defaults + [g_sel_start:120, genome_set_id: 10006],
    // sp_o01: sp_defaults + [num_origins: 1, genome_set_id: 10007],
    // sp_o03: sp_defaults + [num_origins: 3, genome_set_id: 10008],
    // sp_o27: sp_defaults + [num_origins: 27, genome_set_id: 10009],

    // there is bug is previous version where sp_rel and sp_rels03 both use s = 0.3
    // sp_rel:    sp_defaults + [sim_relatedness: 1, s: 0.3, genome_set_id: 30000],
    // sp_rels00: sp_defaults + [sim_relatedness: 1, s: 0.0, genome_set_id: 30002],

    // FOR MODEL exploration

    // change inbreeding starting time (inbreeding transform parameter - power 10: more weights on closer relatives)
    sp_rels03_11a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 10.0, sim_relatedness_g: 200, genome_set_id: 111000],
    sp_rels00_11a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 10.0, sim_relatedness_g: 200, genome_set_id: 110001],
    sp_rels03_11b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 10.0, sim_relatedness_g: 100, genome_set_id: 110002],
    sp_rels00_11b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 10.0, sim_relatedness_g: 100, genome_set_id: 110003],
    sp_rels03_11c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 10.0, sim_relatedness_g:  50, genome_set_id: 110004],
    sp_rels00_11c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 10.0, sim_relatedness_g:  50, genome_set_id: 110005], /////////////
    sp_rels03_12a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 1.0,  sim_relatedness_g: 200, genome_set_id: 120000],
    sp_rels00_12a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 1.0,  sim_relatedness_g: 200, genome_set_id: 120001],
    sp_rels03_12b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 1.0,  sim_relatedness_g: 100, genome_set_id: 120002],
    sp_rels00_12b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 1.0,  sim_relatedness_g: 100, genome_set_id: 120003],
    sp_rels03_12c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 1.0,  sim_relatedness_g:  50, genome_set_id: 120004],
    sp_rels00_12c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 1.0,  sim_relatedness_g:  50, genome_set_id: 120005], ////////////
    sp_rels03_13a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 130000],
    sp_rels00_13a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 130001],
    sp_rels03_13b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 130002],
    sp_rels00_13b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 130003],
    sp_rels03_13c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 130004],
    sp_rels00_13c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 130005],//////////////
    sp_rels03_14a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g: 200, genome_set_id: 140000],
    sp_rels00_14a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g: 200, genome_set_id: 140001],
    sp_rels03_14b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g: 100, genome_set_id: 140002],
    sp_rels00_14b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g: 100, genome_set_id: 140003],
    sp_rels03_14c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g:  50, genome_set_id: 140004],
    sp_rels00_14c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g:  50, genome_set_id: 140005],//////////////

    // change inbreeding starting time (inbreeding transform parameter - power 1: linear)
    sp_rels03_21a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 10.0, sim_relatedness_g: 200, genome_set_id: 211000],
    sp_rels00_21a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 10.0, sim_relatedness_g: 200, genome_set_id: 210001],
    sp_rels03_21b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 10.0, sim_relatedness_g: 100, genome_set_id: 210002],
    sp_rels00_21b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 10.0, sim_relatedness_g: 100, genome_set_id: 210003],
    sp_rels03_21c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 10.0, sim_relatedness_g:  50, genome_set_id: 210004],
    sp_rels00_21c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 10.0, sim_relatedness_g:  50, genome_set_id: 210005], /////////////
    sp_rels03_22a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 1.0,  sim_relatedness_g: 200, genome_set_id: 220000],
    sp_rels00_22a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 1.0,  sim_relatedness_g: 200, genome_set_id: 220001],
    sp_rels03_22b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 1.0,  sim_relatedness_g: 100, genome_set_id: 220002],
    sp_rels00_22b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 1.0,  sim_relatedness_g: 100, genome_set_id: 220003],
    sp_rels03_22c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 1.0,  sim_relatedness_g:  50, genome_set_id: 220004],
    sp_rels00_22c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 1.0,  sim_relatedness_g:  50, genome_set_id: 220005], ////////////
    sp_rels03_23a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 230000],
    sp_rels00_23a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 230001],
    sp_rels03_23b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 230002],
    sp_rels00_23b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 230003],
    sp_rels03_23c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 230004],
    sp_rels00_23c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 230005],//////////////
    sp_rels03_24a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 1e-8, sim_relatedness_g: 200, genome_set_id: 240000],
    sp_rels00_24a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 1e-8, sim_relatedness_g: 200, genome_set_id: 240001],
    sp_rels03_24b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 1e-8, sim_relatedness_g: 100, genome_set_id: 240002],
    sp_rels00_24b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 1e-8, sim_relatedness_g: 100, genome_set_id: 240003],
    sp_rels03_24c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1.0, sim_relatedness_delta: 1e-8, sim_relatedness_g:  50, genome_set_id: 240004],
    sp_rels00_24c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1.0, sim_relatedness_delta: 1e-8, sim_relatedness_g:  50, genome_set_id: 240005],//////////////

    // change inbreeding starting time (inbreeding transform parameter - power 0.1: more weights on relatives, but will put more weights on very close relatives)
    sp_rels03_31a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 10.0, sim_relatedness_g: 200, genome_set_id: 311000],
    sp_rels00_31a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 10.0, sim_relatedness_g: 200, genome_set_id: 310001],
    sp_rels03_31b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 10.0, sim_relatedness_g: 100, genome_set_id: 310002],
    sp_rels00_31b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 10.0, sim_relatedness_g: 100, genome_set_id: 310003],
    sp_rels03_31c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 10.0, sim_relatedness_g:  50, genome_set_id: 310004],
    sp_rels00_31c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 10.0, sim_relatedness_g:  50, genome_set_id: 310005], /////////////
    sp_rels03_32a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 1.0,  sim_relatedness_g: 200, genome_set_id: 320000],
    sp_rels00_32a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 1.0,  sim_relatedness_g: 200, genome_set_id: 320001],
    sp_rels03_32b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 1.0,  sim_relatedness_g: 100, genome_set_id: 320002],
    sp_rels00_32b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 1.0,  sim_relatedness_g: 100, genome_set_id: 320003],
    sp_rels03_32c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 1.0,  sim_relatedness_g:  50, genome_set_id: 320004],
    sp_rels00_32c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 1.0,  sim_relatedness_g:  50, genome_set_id: 320005], ////////////
    sp_rels03_33a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 330000],
    sp_rels00_33a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 330001],
    sp_rels03_33b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 330002],
    sp_rels00_33b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 330003],
    sp_rels03_33c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 330004],
    sp_rels00_33c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 330005],//////////////
    sp_rels03_34a: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 1e-8, sim_relatedness_g: 200, genome_set_id: 340000],
    sp_rels00_34a: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 1e-8, sim_relatedness_g: 200, genome_set_id: 340001],
    sp_rels03_34b: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 1e-8, sim_relatedness_g: 100, genome_set_id: 340002],
    sp_rels00_34b: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 1e-8, sim_relatedness_g: 100, genome_set_id: 340003],
    sp_rels03_34c: sp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 0.1, sim_relatedness_delta: 1e-8, sim_relatedness_g:  50, genome_set_id: 340004],
    sp_rels00_34c: sp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 0.1, sim_relatedness_delta: 1e-8, sim_relatedness_g:  50, genome_set_id: 340005],//////////////

]

def mp_sets = [
    mp_s00: mp_defaults + [s:0.0, genome_set_id: 20000],
    // mp_s01: mp_defaults + [s:0.1, genome_set_id: 20001],
    // mp_s02: mp_defaults + [s:0.2, genome_set_id: 20002],
    mp_s03: mp_defaults + [s:0.3, genome_set_id: 20003],


    // there is bug is previous version where mp_rel and mp_rels03 both use s = 0.3
    // mp_rel:    mp_defaults + [sim_relatedness: 1, s: 0.3, genome_set_id: 30001],
    // mp_rels00: mp_defaults + [sim_relatedness: 1, s: 0.0, genome_set_id: 30003],

    // FOR MODEL exploration
    mp_rels03a: mp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 70000],
    mp_rels00a: mp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 70001],
    mp_rels03b: mp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 70002],
    mp_rels00b: mp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 70003],
    mp_rels03c: mp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 1, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 70004],
    mp_rels00c: mp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 1, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 70005],

    // change inbreeding starting time (inbreeding transform parameter - power 10: more weights on closer relatives)
    mp_rels03d: mp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 80000],
    mp_rels00d: mp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g: 200, genome_set_id: 80001],
    mp_rels03e: mp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 80002],
    mp_rels00e: mp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g: 100, genome_set_id: 80003],
    mp_rels03f: mp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 80004],
    mp_rels00f: mp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 0.01, sim_relatedness_g:  50, genome_set_id: 80005],

    // change inbreeding starting time (inbreeding transform parameter - power 10: more weights on closer relatives, inbreeding transform parameter -- delta: smaller)
    mp_rels03g: mp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g: 200, genome_set_id: 90000],
    mp_rels00g: mp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g: 200, genome_set_id: 90001],
    mp_rels03h: mp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g: 100, genome_set_id: 90002],
    mp_rels00h: mp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g: 100, genome_set_id: 90003],
    mp_rels03i: mp_defaults + [sim_relatedness: 1, s: 0.3, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g:  50, genome_set_id: 90004],
    mp_rels00i: mp_defaults + [sim_relatedness: 1, s: 0.0, sim_relatedness_power: 10, sim_relatedness_delta: 1e-8, sim_relatedness_g:  50, genome_set_id: 90005],
]


process SIM_SP_CHR {
    tag "${args.genome_set_id}_${chrno}"

    publishDir "${resdir}/${args.genome_set_id}_${label}/trees/", pattern: "*.trees", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/vcf/", pattern: "*.vcf.gz", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/daf/", pattern: "*.daf", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/restart_count/", pattern: "*.restart_count", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/true_ne/", pattern: "*.true_ne", mode: 'symlink'

    input: 
    tuple val(label), val(chrno), val(args)

    output: 
    tuple val(label), val(chrno), path("*.trees"), path("*.vcf.gz"), emit: trees_vcf
    tuple val(label), val(chrno), path("*.daf"), emit: daf
    tuple val(label), val(chrno), path("*.restart_count"), emit: restart_count
    tuple val(label), val(chrno), path("*.true_ne"), emit: true_ne

    script:
    def args_chr = (args + [chrno: chrno]).collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    sim_single_pop.py $args_chr

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix="${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,.true_ne}
    """
}

// for inbreeding simulation
process SIM_SP_GW {
    tag "${args.genome_set_id}_genome_wide"

    publishDir "${resdir}/${args.genome_set_id}_${label}/trees/", pattern: "*.trees", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/vcf/", pattern: "*.vcf.gz", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/daf/", pattern: "*.daf", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/restart_count/", pattern: "*.restart_count", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/true_ne/", pattern: "*.true_ne", mode: 'symlink'

    input: 
    tuple val(label), val(args)

    output: 
    tuple val(label), path("*.trees"), path("*.vcf.gz"), emit: trees_vcf
    tuple val(label), path("*.daf"), emit: daf
    tuple val(label), path("*.restart_count"), emit: restart_count
    tuple val(label), path("*.true_ne"), emit: true_ne

    script:
    def args_chr = args.collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    sim_single_pop_whole_genome.py $args_chr

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix="${args.genome_set_id}"
    """
    touch ${prefix}_{1..14}{.trees,.vcf.gz,.daf,.restart_count,.true_ne}
    """
}

process SIM_MP_CHR {
    tag "${args.genome_set_id}_${chrno}"

    publishDir "${resdir}/${args.genome_set_id}_${label}/trees/", pattern: "*.trees", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/vcf/", pattern: "*.vcf.gz", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/restart_count/", pattern: "*.restart_count", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/daf/", pattern: "*.daf", mode: 'symlink'
    publishDir "${resdir}/${args.genome_set_id}_${label}/demog/", pattern: "*_demog.png", mode: 'symlink'

    input:
    tuple val(label), val(chrno), val(args)

    output:
    tuple val(label), val(chrno), path("*.trees"), path("*.vcf.gz"), emit: trees_vcf
    tuple val(label), val(chrno), path("*.daf"), emit: daf
    tuple val(label), val(chrno), path("*.restart_count"), emit: restart_count
    tuple val(label), val(chrno), path("*_demog.png"), emit: demog

    script:
    def args_chr = (args + [chrno: chrno]).collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    sim_multiple_pop.py $args_chr

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix="${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,_demog.png}
    """
}

process CALL_IBD {
    tag "${args.genome_set_id}_${chrno}"

    publishDir "${resdir}/${args.genome_set_id}_${label}/tskibd/", pattern: "*_tskibd.ibd", mode: 'symlink'

    input: 
        tuple val(label), val(chrno), val(args),  path(trees), path(vcf)
    output:
        tuple val(label), val(chrno), path("*_tskibd.ibd"), emit: tskibd
    script:
    def args_local = [
        tree: trees,
        vcf: vcf,
        chrno: chrno,
        r: args.r,
        seqlen: args.seqlen,
        genome_set_id: args.genome_set_id,
        mincm: 2.0, call_hmmibd: 0, n:100, m:5,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    call_ibd.py ${args_local}
    """
    stub:
    def prefix="${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.map,_tskibd.ibd}
    """
}

process PROC_DIST_NE {
    tag "${genome_set_id}"

    publishDir "${resdir}/${genome_set_id}_${label}/ne_input/", pattern: "*.sh", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/ne_input/", pattern: "*.map", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/ne_input/", pattern: "*.ibd.gz", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/ibddist_ibd/", pattern: "*.ibddist.ibdobj.gz", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/ibdne_ibd/", pattern: "*.ibdne.ibdobj.gz", mode: 'symlink'

    input:
        tuple val(label), path(ibd_lst), path(vcf_lst), val(genome_set_id)
    output:
        tuple val(label), path("ibdne.jar"), path("*_orig.sh"), \
                path("*_orig.map"), path("*_orig.ibd.gz"), emit: ne_input_orig
        tuple val(label), path("ibdne.jar"), path("*_rmpeaks.sh"),  \
                path("*_rmpeaks.map"), path("*_rmpeaks.ibd.gz"), emit: ne_input_rmpeaks
        tuple val(label), path("*.ibddist.ibdobj.gz"), emit: ibddist_ibd_obj
        tuple val(label), path("*.ibdne.ibdobj.gz"), emit: ibdne_ibd_obj
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        vcf_files: "${vcf_lst}", // path is a blank separate list
        genome_set_id: genome_set_id,
        ibdne_mincm: params.ibdne_mincm,
        ibdne_minregion: params.ibdne_minregion,
        ibdne_flatmeth: params.ibdne_flatmeth,
        ibdne_no_diploid_conversion: params.ibdne_no_diploid_convertion,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")

    """
    proc_dist_ne.py ${args_local} 
    """
    stub:
    """
    touch ibdne.jar
    touch ${genome_set_id}{_orig.sh,_orig.map,_orig.ibd.gz}
    touch ${genome_set_id}{_rmpeaks.sh,_rmpeaks.map,_rmpeaks.ibd.gz}
    touch ${genome_set_id}.ibddist.ibdobj.gz
    touch ${genome_set_id}_orig.ibdne.ibdobj.gz
    touch ${genome_set_id}_rmpeaks.ibdne.ibdobj.gz
    """
}

process PROC_INFOMAP {
    tag "${genome_set_id}"

    publishDir "${resdir}/${genome_set_id}_${label}/ifm_input/", \
        pattern: "*.ibdobj.gz", mode: 'symlink'

    input:
        tuple val(label), path(ibd_lst), path(vcf_lst), val(genome_set_id)
    output:
        tuple val(label), path("*_orig.ifm.ibdobj.gz"), emit: ifm_orig_ibd_obj
        tuple val(label), path("*_rmpeaks.ifm.ibdobj.gz"), emit: ifm_rmpeaks_ibd_obj
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        vcf_files: "${vcf_lst}", // path is a blank separate list
        genome_set_id: genome_set_id,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    proc_infomap.py ${args_local}
    """
    stub:
    """
    touch ${genome_set_id}{_orig.ifm.ibdobj.gz,_rmpeaks.ifm.ibdobj.gz}
    """
}

process RUN_IBDNE {
    tag "${args.genome_set_id}_${are_peaks_removed}"

    publishDir "${resdir}/${args.genome_set_id}_${label}/ne_output/",  mode: 'symlink'

    input:
        tuple val(label), path(ibdne_jar), path(ibdne_sh), path(gmap), path(ibd_gz), \
            val(are_peaks_removed), val(args)
    output:
        tuple val(label), val(are_peaks_removed), path("*.ne")
    script:
    """
    bash ${ibdne_sh}
    """
    stub:
    def src = are_peaks_removed ? "rmpeaks": "orig"
    """
    touch ${args.genome_set_id}_${src}.ne
    """
}

process RUN_INFOMAP {
    tag "${args.genome_set_id}_${are_peaks_removed}"
    publishDir "${resdir}/${args.genome_set_id}_${label}/ifm_output/",  mode: 'symlink'
    input:
        tuple val(label), path(ibd_obj), val(are_peaks_removed), val(args)
    output:
        tuple val(label), val(are_peaks_removed), path("*_member.pq")
    script:
    def cut_mode = are_peaks_removed? 'rmpeaks': 'orig'
    def args_local = [
        ibd_obj: ibd_obj,
        npop: args.npop,
        nsam: args.nsam,
        genome_set_id: args.genome_set_id,
        cut_mode: cut_mode,
        ntrials: params.ifm_ntrials,
        transform: params.ifm_transform,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    run_infomap.py ${args_local}
    """
    stub:
    def cut_mode = are_peaks_removed? 'rmpeaks': 'orig'
    """
    touch ${args.genome_set_id}_${cut_mode}_member.pq
    """
}



workflow WF_SP {
    // SIMULATION SETS

    // println("\n\n------------- Single Population Sets-------------")
    // sp_sets.each{k, v -> println("${k}\t\t${v}")}

    def expanded_sets = [:]

    if(params.num_reps == 1) {
        expanded_sets = sp_sets
    } else{
        (0..<(params.num_reps)).each {rep ->
            sp_sets.each{k, v->
                def k1 = "${k}_rep${rep}"
                def v1 = v + [genome_set_id: v.genome_set_id + 10 *  rep]
                expanded_sets[k1] = v1
            }
        } 
    }



    // SIMULATE chromosomes
    chr_chrno = channel.from(1..14)
    ch_sp_params = channel.from(expanded_sets.collect{k, v-> [k, v]})

    if (params.test) {
        ch_sp_params = ch_sp_params.first().map{
            label, args -> def args2 = args + [nsam:50]; [label, args2]
        }
    }

    if (params.sim_inbreeding) {
        SIM_SP_GW(ch_sp_params)
        // need sort by chromosome numbers
        ch_trees = SIM_SP_GW.out.trees_vcf.flatMap{label, trees, vcf -> 
            trees.collect{t-> [label, t.getSimpleName().split("_")[-1].toInteger(), t]}
        }
        ch_vcf = SIM_SP_GW.out.trees_vcf.flatMap{label, trees, vcf ->
            vcf.collect{v-> [label, v.getSimpleName().split("_")[-1].toInteger(), v]}
        } 
        ch_trees_vcf = ch_trees.combine(ch_vcf, by:[0,1])
    } else{
        SIM_SP_CHR(ch_sp_params.combine(chr_chrno).map{a,b,c->[a,c,b]})
        ch_trees_vcf = SIM_SP_CHR.out.trees_vcf
    }
    
    
//    ch_trees_vcf.combine(ch_sp_params, by:0)
 //       .map{label,chrno,trees,vcf,args-> [label,chrno, args, trees, vcf]}.take(1).view()

    // CALL IBD segment per chromosome
    CALL_IBD(ch_trees_vcf.combine(ch_sp_params, by:0)
        .map{label,chrno,trees,vcf,args-> [label,chrno, args, trees, vcf]}
    )
    // CALL_IBD.out.tskibd.view{label, chrno, ibd -> [label, chrno, ibd.getName()]}



    // COLLECT PER SETS
    ch_ibd_per_genome = CALL_IBD.out.tskibd
        // tskibd: label, chrno, ibd
        // trees_vcf: label, chrno, tree, vcf
        .combine(ch_trees_vcf.map{it[[0, 1, 3]]}, by: [0, 1])
        .map{ label, chrno, ibd, vcf -> 
            [groupKey(label, 14), [chrno, ibd, vcf]] }
        .groupTuple(by:0, sort: {x,y-> x[0]<=>y[0] })
        .combine(ch_sp_params, by:0)
        .map{label, ll, args-> 
            [label, ll.collect{it[1]}, ll.collect{it[2]}, args.genome_set_id]}

    // ch_ibd_per_genome.view{label, trees_lst, genome_set_id->
    //        [label, trees_lst.collect{it.getName()}, genome_set_id]}


    // Process IBD for ibd distribution and ne analyses
    PROC_DIST_NE(ch_ibd_per_genome)
    // PROC_DIST_NE.out.ibddist_ibd_obj.view{label, ibdpq -> [label, ibdpq.getName()]}


    // RUN IBDNe actually
    RUN_IBDNE(
        PROC_DIST_NE.out.ne_input_orig.map{it -> it + false}.map{it-> it + expanded_sets[it[0]]}
        // use mix so the second part do not have to wait for first part to finish complely
        .mix(
            PROC_DIST_NE.out.ne_input_rmpeaks.map{it -> it + true}.map{it-> it + expanded_sets[it[0]]}
        )
    )

    // RUN_IBDNE.out.map{label, yn, ne->[label, yn, ne.getName()]}.view()
    

}


workflow WF_MP {
    // SIMULATION SETS
    
    // println("\n\n------------- Multiple Population Sets-------------")
    // mp_sets.each{k, v -> println("${k}\t\t${v}")}

    def expanded_sets = [:]

    if(params.num_reps == 1)
    {
        expanded_sets = mp_sets
    }
    else {
        (0..<(params.num_reps)).each {rep ->
            mp_sets.each{k, v->
                def k1 = "${k}_rep${rep}"
                def v1 = v + [genome_set_id: v.genome_set_id + 10 *  rep]
                expanded_sets[k1] = v1
            }
        } 

    }

    // SIMULATE chromosomes
    chr_chrno = channel.from(1..14)
    ch_mp_params = channel.from(expanded_sets.collect{k, v-> [k, v]})

    if (params.test) {
        ch_mp_params = ch_mp_params.first().map{
            label, args -> def args2 = args + [nsam:50, npop:2, N:500]; [label, args2]
        }
    }

    SIM_MP_CHR(ch_mp_params.combine(chr_chrno).map{a,b,c->[a,c,b]})
    // SIM_MP_CHR.out.trees_vcf.view{it-> 
    //    "${it[0]}\t${it[1]}\t${it[2].getName()}\t${it[3].getName()}"}



    // CALL IBD segment per chromosome
    CALL_IBD(SIM_MP_CHR.out.trees_vcf.combine(ch_mp_params, by:0)
        .map{label,chrno,trees,vcf,args-> [label,chrno, args, trees, vcf]}
    )
    // CALL_IBD.out.tskibd.view{label, chrno, ibd -> [label, chrno, ibd.getName()]}


    // COLLECT PER SETS
    ch_ibd_per_genome = CALL_IBD.out.tskibd
        // tskibd: label, chrno, ibd
        // trees_vcf: label, chrno, tree, vcf
        .combine(SIM_MP_CHR.out.trees_vcf.map{it[[0, 1, 3]]}, by: [0, 1])
        .map{ label, chrno, ibd, vcf -> 
            [groupKey(label, 14), [chrno, ibd, vcf]] }
        .groupTuple(by:0, sort: {x,y-> x[0]<=>y[0] })
        .combine(ch_mp_params, by:0)
        .map{label, ll, args-> 
            [label, ll.collect{it[1]}, ll.collect{it[2]}, args.genome_set_id]}

    // ch_ibd_per_genome.view{label, trees_lst, genome_set_id->
    //        [label, trees_lst.collect{it.getName()}, genome_set_id]}


    // Process IBD for ibd distribution and ne analyses
    PROC_INFOMAP(ch_ibd_per_genome)
    
    // PROC_INFOMAP.out.ifm_orig_ibd_ob.view{label, ibdpq -> [label, ibdpq.getName()]}


    // RUN INFOMAP
    RUN_INFOMAP(
        PROC_INFOMAP.out.ifm_orig_ibd_obj.map{it -> it + false}.combine(ch_mp_params, by:0).concat(
            PROC_INFOMAP.out.ifm_rmpeaks_ibd_obj.map{it -> it + true}.combine(ch_mp_params, by:0)
        )
    )

    // RUN_INFOMAP.out.view{label, yns, member -> [label, yns, member.getName()]}

}

workflow {
    WF_SP()
    WF_MP()
}
