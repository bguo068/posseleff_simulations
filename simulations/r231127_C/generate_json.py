import json
from pathlib import Path
from typing import Any

sp_shared: dict[str, Any] = dict(
    N=10000,
    N0=1000,
    g_sel_start=80,
    sim_relatedness_g=200,
    sim_relatedness_power=16,
    sim_relatedness_bypass=4,
    sim_relatedness_bypass_complement=1,
)


def main():
    res = {}
    sp_label = "sp"
    idx = 10000

    for s in [0.0, 0.3]:
        s_label = f"s{s}"
        for sim_relatedness_delta in [1000, 1.0, 0.3, 0.1, 0.03, 0.01, 0.001]:
            sim_relatedness_delta_label = f"srd{sim_relatedness_delta}"
            if sim_relatedness_delta == 1000:
                label = f"{sp_label}_{s_label}_{sim_relatedness_delta_label}"
                val = sp_shared.copy()
                val.update(
                    {
                        "genome_set_id": idx,
                        "s": s,
                        "sim_relatedness": 0,
                        "sim_relatedness_delta": sim_relatedness_delta,
                    }
                )
                res[label] = val
                idx += 10000
            else:
                label = "_".join(
                    [
                        sp_label,
                        s_label,
                        sim_relatedness_delta_label,
                    ]
                )
                val = sp_shared.copy()
                val.update(
                    {
                        "genome_set_id": idx,
                        "s": s,
                        "sim_relatedness": 1,
                        "sim_relatedness_delta": sim_relatedness_delta,
                    }
                )
                res[label] = val
                idx += 10000

    json_str = json.dumps(res, indent=4)
    Path("sp_genome_sets.json").write_text(json_str)


if __name__ == "__main__":
    main()
