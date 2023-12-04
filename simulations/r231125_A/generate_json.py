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
                for N0 in [10000, 3000, 1000, 300]:
                    if N0 == 10000:
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
                        n0_label = f"n0{N0}"
                        for g_ne_change_start in [10, 25, 50, 100, 200]:
                            g_ne_change_start_label = f"ssg{g_ne_change_start}"
                            label = "_".join(
                                [
                                    mp_label,
                                    sel_mig_label,
                                    s_label,
                                    g_sel_start_label,
                                    n0_label,
                                    g_ne_change_start_label,
                                ]
                            )
                            val = mp_shared.copy()
                            val.update(
                                {
                                    "genome_set_id": idx,
                                    "sel_mig": sel_mig,
                                    "s": s,
                                    "g_sel_start": g_sel_start,
                                    "N0": N0,
                                    "g_ne_change_start": g_ne_change_start,
                                }
                            )
                            res[label] = val
                            idx += 10000

        json_str = json.dumps(res, indent=4)
        Path("mp_genome_sets.json").write_text(json_str)


if __name__ == "__main__":
    main()
