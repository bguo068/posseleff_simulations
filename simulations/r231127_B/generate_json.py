from pathlib import Path
import json
from typing import Any

sp_shared: dict[str, Any] = dict(N=10000)


def main():
    res = {}
    sp_label = "sp"
    idx = 100
    sim_selfing_g = 200

    for s in [0.0, 0.3]:
        for sim_selfing_rate in [0.0, 0.1, 0.2, 0.3, 0.5, 0.7]:
            val = dict(
                s=s,
                sim_selfing_rate=sim_selfing_rate,
                genome_set_id=idx,
                sim_selfing_g=sim_selfing_g,
            )
            label = f"{sp_label}_s{s}_ssr{sim_selfing_rate}"
            res[label] = val
            # print(label, val)

            idx += 100
    # print(res)

    json_str = json.dumps(res, indent=4)

    Path("sp_genome_sets.json").write_text(json_str)


main()
