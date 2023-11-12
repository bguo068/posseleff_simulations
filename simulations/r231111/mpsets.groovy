import groovy.json.JsonOutput

static void main(String[] args2) {
    def res = [:]
    def mp_label = "mp"
    def idx = 10000
    for (s in [0.0, 0.3]) {
        def s_label = s == 0.3 ? "s03" : "s00"
        for (g_sel_start in [50, 80]) {
            def g_sel_start_label = "gss${g_sel_start}"
            for (sim_relatedness in [0, 1]) {
                if (sim_relatedness == 0) {
                    def label = [mp_label, s_label, g_sel_start_label].join("_")
                    res[label] = [genome_set_id: idx, s: s, g_sel_start: g_sel_start, sim_relatedness: 0]
                    idx += 10000
                } else {
                    for (sim_relatedness_power in [16, 8]) {
                        def sim_relatedness_power_label = "amp${sim_relatedness_power}"
                        for (sim_relatedness_delta in [0.5, 0.3, 0.1, 0.01, 0.001]) {
                            // def d = [0.01: 2, 0.001: 3]
                            def d = [0.01: 2, 0.001: 3, 0.5: 4, 0.3: 5, 0.1: 6]
                            def sim_relatedness_delta_label = "amd${d[sim_relatedness_delta]}"
                            for (sim_relatedness_g in [100, 50, 25]) {
                                def sim_relatedness_g_label = "amg${sim_relatedness_g}"
                                def label = [mp_label, s_label, g_sel_start_label, sim_relatedness_power_label, sim_relatedness_delta_label, sim_relatedness_g_label].join("_")
                                res[label] = [genome_set_id        : idx,
                                              s                    : s, g_sel_start: g_sel_start,
                                              sim_relatedness      : 1,
                                              sim_relatedness_power: sim_relatedness_power,
                                              sim_relatedness_delta: sim_relatedness_delta,
                                              sim_relatedness_g    : sim_relatedness_g,]
                                idx += 10000
                            }
                        }
                    }
                }
            }
        }
    }
//    print(res)
    /*
    def x = [
            mp_s00_gss50                  : [genome_set_id: 10000, s: 0.0, g_sel_start: 50, sim_relatedness: 0],
            mp_s00_gss50_amp16_amd4_amg100: [genome_set_id: 20000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 100],
            mp_s00_gss50_amp16_amd4_amg50 : [genome_set_id: 30000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 50],
            mp_s00_gss50_amp16_amd4_amg25 : [genome_set_id: 40000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 25],
            mp_s00_gss50_amp16_amd5_amg100: [genome_set_id: 50000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 100],
            mp_s00_gss50_amp16_amd5_amg50 : [genome_set_id: 60000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 50],
            mp_s00_gss50_amp16_amd5_amg25 : [genome_set_id: 70000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 25],
            mp_s00_gss50_amp16_amd6_amg100: [genome_set_id: 80000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 100],
            mp_s00_gss50_amp16_amd6_amg50 : [genome_set_id: 90000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 50],
            mp_s00_gss50_amp16_amd6_amg25 : [genome_set_id: 100000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 25],
            mp_s00_gss50_amp16_amd2_amg100: [genome_set_id: 110000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 100],
            mp_s00_gss50_amp16_amd2_amg50 : [genome_set_id: 120000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 50],
            mp_s00_gss50_amp16_amd2_amg25 : [genome_set_id: 130000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 25],
            mp_s00_gss50_amp16_amd3_amg100: [genome_set_id: 140000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 100],
            mp_s00_gss50_amp16_amd3_amg50 : [genome_set_id: 150000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 50],
            mp_s00_gss50_amp16_amd3_amg25 : [genome_set_id: 160000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 25],
            mp_s00_gss50_amp8_amd4_amg100 : [genome_set_id: 170000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 100],
            mp_s00_gss50_amp8_amd4_amg50  : [genome_set_id: 180000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 50],
            mp_s00_gss50_amp8_amd4_amg25  : [genome_set_id: 190000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 25],
            mp_s00_gss50_amp8_amd5_amg100 : [genome_set_id: 200000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 100],
            mp_s00_gss50_amp8_amd5_amg50  : [genome_set_id: 210000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 50],
            mp_s00_gss50_amp8_amd5_amg25  : [genome_set_id: 220000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 25],
            mp_s00_gss50_amp8_amd6_amg100 : [genome_set_id: 230000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 100],
            mp_s00_gss50_amp8_amd6_amg50  : [genome_set_id: 240000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 50],
            mp_s00_gss50_amp8_amd6_amg25  : [genome_set_id: 250000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 25],
            mp_s00_gss50_amp8_amd2_amg100 : [genome_set_id: 260000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 100],
            mp_s00_gss50_amp8_amd2_amg50  : [genome_set_id: 270000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 50],
            mp_s00_gss50_amp8_amd2_amg25  : [genome_set_id: 280000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 25],
            mp_s00_gss50_amp8_amd3_amg100 : [genome_set_id: 290000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 100],
            mp_s00_gss50_amp8_amd3_amg50  : [genome_set_id: 300000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 50],
            mp_s00_gss50_amp8_amd3_amg25  : [genome_set_id: 310000, s: 0.0, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 25],
            mp_s00_gss80                  : [genome_set_id: 320000, s: 0.0, g_sel_start: 80, sim_relatedness: 0],
            mp_s00_gss80_amp16_amd4_amg100: [genome_set_id: 330000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 100],
            mp_s00_gss80_amp16_amd4_amg50 : [genome_set_id: 340000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 50],
            mp_s00_gss80_amp16_amd4_amg25 : [genome_set_id: 350000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 25],
            mp_s00_gss80_amp16_amd5_amg100: [genome_set_id: 360000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 100],
            mp_s00_gss80_amp16_amd5_amg50 : [genome_set_id: 370000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 50],
            mp_s00_gss80_amp16_amd5_amg25 : [genome_set_id: 380000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 25],
            mp_s00_gss80_amp16_amd6_amg100: [genome_set_id: 390000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 100],
            mp_s00_gss80_amp16_amd6_amg50 : [genome_set_id: 400000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 50],
            mp_s00_gss80_amp16_amd6_amg25 : [genome_set_id: 410000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 25],
            mp_s00_gss80_amp16_amd2_amg100: [genome_set_id: 420000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 100],
            mp_s00_gss80_amp16_amd2_amg50 : [genome_set_id: 430000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 50],
            mp_s00_gss80_amp16_amd2_amg25 : [genome_set_id: 440000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 25],
            mp_s00_gss80_amp16_amd3_amg100: [genome_set_id: 450000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 100],
            mp_s00_gss80_amp16_amd3_amg50 : [genome_set_id: 460000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 50],
            mp_s00_gss80_amp16_amd3_amg25 : [genome_set_id: 470000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 25],
            mp_s00_gss80_amp8_amd4_amg100 : [genome_set_id: 480000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 100],
            mp_s00_gss80_amp8_amd4_amg50  : [genome_set_id: 490000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 50],
            mp_s00_gss80_amp8_amd4_amg25  : [genome_set_id: 500000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 25],
            mp_s00_gss80_amp8_amd5_amg100 : [genome_set_id: 510000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 100],
            mp_s00_gss80_amp8_amd5_amg50  : [genome_set_id: 520000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 50],
            mp_s00_gss80_amp8_amd5_amg25  : [genome_set_id: 530000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 25],
            mp_s00_gss80_amp8_amd6_amg100 : [genome_set_id: 540000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 100],
            mp_s00_gss80_amp8_amd6_amg50  : [genome_set_id: 550000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 50],
            mp_s00_gss80_amp8_amd6_amg25  : [genome_set_id: 560000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 25],
            mp_s00_gss80_amp8_amd2_amg100 : [genome_set_id: 570000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 100],
            mp_s00_gss80_amp8_amd2_amg50  : [genome_set_id: 580000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 50],
            mp_s00_gss80_amp8_amd2_amg25  : [genome_set_id: 590000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 25],
            mp_s00_gss80_amp8_amd3_amg100 : [genome_set_id: 600000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 100],
            mp_s00_gss80_amp8_amd3_amg50  : [genome_set_id: 610000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 50],
            mp_s00_gss80_amp8_amd3_amg25  : [genome_set_id: 620000, s: 0.0, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 25],
            mp_s03_gss50                  : [genome_set_id: 630000, s: 0.3, g_sel_start: 50, sim_relatedness: 0],
            mp_s03_gss50_amp16_amd4_amg100: [genome_set_id: 640000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 100],
            mp_s03_gss50_amp16_amd4_amg50 : [genome_set_id: 650000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 50],
            mp_s03_gss50_amp16_amd4_amg25 : [genome_set_id: 660000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 25],
            mp_s03_gss50_amp16_amd5_amg100: [genome_set_id: 670000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 100],
            mp_s03_gss50_amp16_amd5_amg50 : [genome_set_id: 680000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 50],
            mp_s03_gss50_amp16_amd5_amg25 : [genome_set_id: 690000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 25],
            mp_s03_gss50_amp16_amd6_amg100: [genome_set_id: 700000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 100],
            mp_s03_gss50_amp16_amd6_amg50 : [genome_set_id: 710000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 50],
            mp_s03_gss50_amp16_amd6_amg25 : [genome_set_id: 720000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 25],
            mp_s03_gss50_amp16_amd2_amg100: [genome_set_id: 730000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 100],
            mp_s03_gss50_amp16_amd2_amg50 : [genome_set_id: 740000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 50],
            mp_s03_gss50_amp16_amd2_amg25 : [genome_set_id: 750000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 25],
            mp_s03_gss50_amp16_amd3_amg100: [genome_set_id: 760000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 100],
            mp_s03_gss50_amp16_amd3_amg50 : [genome_set_id: 770000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 50],
            mp_s03_gss50_amp16_amd3_amg25 : [genome_set_id: 780000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 25],
            mp_s03_gss50_amp8_amd4_amg100 : [genome_set_id: 790000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 100],
            mp_s03_gss50_amp8_amd4_amg50  : [genome_set_id: 800000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 50],
            mp_s03_gss50_amp8_amd4_amg25  : [genome_set_id: 810000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 25],
            mp_s03_gss50_amp8_amd5_amg100 : [genome_set_id: 820000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 100],
            mp_s03_gss50_amp8_amd5_amg50  : [genome_set_id: 830000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 50],
            mp_s03_gss50_amp8_amd5_amg25  : [genome_set_id: 840000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 25],
            mp_s03_gss50_amp8_amd6_amg100 : [genome_set_id: 850000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 100],
            mp_s03_gss50_amp8_amd6_amg50  : [genome_set_id: 860000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 50],
            mp_s03_gss50_amp8_amd6_amg25  : [genome_set_id: 870000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 25],
            mp_s03_gss50_amp8_amd2_amg100 : [genome_set_id: 880000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 100],
            mp_s03_gss50_amp8_amd2_amg50  : [genome_set_id: 890000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 50],
            mp_s03_gss50_amp8_amd2_amg25  : [genome_set_id: 900000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 25],
            mp_s03_gss50_amp8_amd3_amg100 : [genome_set_id: 910000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 100],
            mp_s03_gss50_amp8_amd3_amg50  : [genome_set_id: 920000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 50],
            mp_s03_gss50_amp8_amd3_amg25  : [genome_set_id: 930000, s: 0.3, g_sel_start: 50, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 25],
            mp_s03_gss80                  : [genome_set_id: 940000, s: 0.3, g_sel_start: 80, sim_relatedness: 0],
            mp_s03_gss80_amp16_amd4_amg100: [genome_set_id: 950000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 100],
            mp_s03_gss80_amp16_amd4_amg50 : [genome_set_id: 960000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 50],
            mp_s03_gss80_amp16_amd4_amg25 : [genome_set_id: 970000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.5, sim_relatedness_g: 25],
            mp_s03_gss80_amp16_amd5_amg100: [genome_set_id: 980000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 100],
            mp_s03_gss80_amp16_amd5_amg50 : [genome_set_id: 990000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 50],
            mp_s03_gss80_amp16_amd5_amg25 : [genome_set_id: 1000000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.3, sim_relatedness_g: 25],
            mp_s03_gss80_amp16_amd6_amg100: [genome_set_id: 1010000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 100],
            mp_s03_gss80_amp16_amd6_amg50 : [genome_set_id: 1020000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 50],
            mp_s03_gss80_amp16_amd6_amg25 : [genome_set_id: 1030000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.1, sim_relatedness_g: 25],
            mp_s03_gss80_amp16_amd2_amg100: [genome_set_id: 1040000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 100],
            mp_s03_gss80_amp16_amd2_amg50 : [genome_set_id: 1050000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 50],
            mp_s03_gss80_amp16_amd2_amg25 : [genome_set_id: 1060000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.01, sim_relatedness_g: 25],
            mp_s03_gss80_amp16_amd3_amg100: [genome_set_id: 1070000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 100],
            mp_s03_gss80_amp16_amd3_amg50 : [genome_set_id: 1080000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 50],
            mp_s03_gss80_amp16_amd3_amg25 : [genome_set_id: 1090000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 16, sim_relatedness_delta: 0.001, sim_relatedness_g: 25],
            mp_s03_gss80_amp8_amd4_amg100 : [genome_set_id: 1100000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 100],
            mp_s03_gss80_amp8_amd4_amg50  : [genome_set_id: 1110000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 50],
            mp_s03_gss80_amp8_amd4_amg25  : [genome_set_id: 1120000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.5, sim_relatedness_g: 25],
            mp_s03_gss80_amp8_amd5_amg100 : [genome_set_id: 1130000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 100],
            mp_s03_gss80_amp8_amd5_amg50  : [genome_set_id: 1140000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 50],
            mp_s03_gss80_amp8_amd5_amg25  : [genome_set_id: 1150000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.3, sim_relatedness_g: 25],
            mp_s03_gss80_amp8_amd6_amg100 : [genome_set_id: 1160000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 100],
            mp_s03_gss80_amp8_amd6_amg50  : [genome_set_id: 1170000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 50],
            mp_s03_gss80_amp8_amd6_amg25  : [genome_set_id: 1180000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.1, sim_relatedness_g: 25],
            mp_s03_gss80_amp8_amd2_amg100 : [genome_set_id: 1190000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 100],
            mp_s03_gss80_amp8_amd2_amg50  : [genome_set_id: 1200000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 50],
            mp_s03_gss80_amp8_amd2_amg25  : [genome_set_id: 1210000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.01, sim_relatedness_g: 25],
            mp_s03_gss80_amp8_amd3_amg100 : [genome_set_id: 1220000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 100],
            mp_s03_gss80_amp8_amd3_amg50  : [genome_set_id: 1230000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 50],
            mp_s03_gss80_amp8_amd3_amg25  : [genome_set_id: 1240000, s: 0.3, g_sel_start: 80, sim_relatedness: 1, sim_relatedness_power: 8, sim_relatedness_delta: 0.001, sim_relatedness_g: 25]
    ]
     */


    print(JsonOutput.prettyPrint(JsonOutput.toJson(res)))

}
