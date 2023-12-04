import json
from pathlib import Path

mp_shared = dict(N=10000)


def main():
    res = {}
    mp_label = "mp"
    idx = 10000

    for sel_mig in [0.01, 0.015, 0.03]:
        sel_mig_label = f"smig{sel_mig}"
        for s in [0.0, 0.3]:
            s_label = "s03" if s == 0.3 else "s00"
            for g_sel_start in [80]:
                g_sel_start_label = f"gss{g_sel_start}"
                for sim_selfing_rate in [0.0, 0.1, 0.2, 0.5]:
                    if sim_selfing_rate == 0.0:
                        label = (
                            f"{mp_label}_{sel_mig_label}_{s_label}_{g_sel_start_label}"
                        )
                        val = mp_shared.copy()
                        val.update(
                            {
                                "genome_set_id": idx,
                                "sel_mig": sel_mig,
                                "s": s,
                                "g_sel_start": g_sel_start,
                            }
                        )
                        res[label] = val
                        idx += 10000
                    else:
                        # for sim_relatedness_power in [16, 8]:
                        sim_selfing_rate_label = f"ssr{sim_selfing_rate}"
                        for sim_selfing_g in [10, 25, 50, 100, 200]:
                            sim_selfing_g_label = f"ssg{sim_selfing_g}"
                            label = "_".join(
                                [
                                    mp_label,
                                    sel_mig_label,
                                    s_label,
                                    g_sel_start_label,
                                    sim_selfing_rate_label,
                                    sim_selfing_g_label,
                                ]
                            )
                            val = mp_shared.copy()
                            val.update(
                                {
                                    "genome_set_id": idx,
                                    "sel_mig": sel_mig,
                                    "s": s,
                                    "g_sel_start": g_sel_start,
                                    "sim_selfing_rate": sim_selfing_rate,
                                    "sim_selfing_g": sim_selfing_g,
                                }
                            )
                            res[label] = val
                            idx += 10000

        json_str = json.dumps(res, indent=4)
        Path("mp_genome_sets.json").write_text(json_str)


if __name__ == "__main__":
    main()
