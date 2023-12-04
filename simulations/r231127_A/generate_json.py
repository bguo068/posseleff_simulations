from pathlib import Path
import json
from typing import Any

sp_shared: dict[str, Any] = dict(N=10000)


def main():
    res = {}
    sp_label = "sp"
    idx = 100

    for s in [0.0, 0.3]:
        for N0 in [500, 1000, 3000, 10000]:
            val = dict(
                s=s,
                N0=N0,
                genome_set_id=idx,
            )
            label = f"{sp_label}_s{s}_N0{N0}"
            res[label] = val
            # print(label, val)

            idx += 100
    # print(res)

    json_str = json.dumps(res, indent=4)

    Path("sp_genome_sets.json").write_text(json_str)


main()
