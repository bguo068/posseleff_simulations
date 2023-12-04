import json
from pathlib import Path
from typing import Any

mp_shared: dict[str, Any] = dict(
    N=10000,
    g_sel_start=80,
    sel_mig=0.01,
    g_ne_change_start=200,
    sim_selfing_g=200,
)


def main():
    res = {}
    mp_label = "mp"

    for i, sim_selfing_rate in enumerate([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9]):
        sim_selfing_rate_label = f"ssr{sim_selfing_rate}"
        for j, s in enumerate([0.0, 0.3]):
            s_label = f"s{s}"
            for k, N0 in enumerate([10000, 1000]):
                if N0 == 10000:
                    label = f"{mp_label}_{s_label}_{sim_selfing_rate_label}"
                    val = mp_shared.copy()
                    val.update(
                        {
                            "genome_set_id": i * 1000000 + j * 100000 + k * 100,
                            "s": s,
                            "sim_selfing_rate": sim_selfing_rate,
                        }
                    )
                    res[label] = val
                else:
                    n0_label = f"n0{N0}"
                    label = "_".join(
                        [
                            mp_label,
                            s_label,
                            n0_label,
                            sim_selfing_rate_label,
                        ]
                    )
                    val = mp_shared.copy()
                    val.update(
                        {
                            "genome_set_id": i * 1000000 + j * 100000 + k * 100,
                            "s": s,
                            "N0": N0,
                            "sim_selfing_rate": sim_selfing_rate,
                        }
                    )
                    res[label] = val

    json_str = json.dumps(res, indent=4)
    Path("mp_genome_sets.json").write_text(json_str)


if __name__ == "__main__":
    main()
