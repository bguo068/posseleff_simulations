#! /usr/bin/env python3
import os
import json
import sys
from pathlib import Path
from subprocess import check_output, run

# Install conda environment
def is_conda_available():
    return Path(os.environ["CONDA_EXE"]).exists()


def is_simulation_env_installed():
    jstr = check_output("conda env list --json", shell=True, text=True)
    envs = json.loads(jstr)["envs"]
    res = False
    for env in envs:
        if Path(env).name == "simulation":
            res = True
            break
    return res


def install_simulation_env(env_file: str):
    run(f"mamba env create -f {env_file}", shell=True, check=True)


# Install tskibd
def install_tskibd():
    curdir = Path(".").absolute()
    (curdir / "bin").mkdir(parents=True, exist_ok=True)
    if not (curdir / "bin/tskibd").exists():
        run(
            f"""
            eval "$(conda shell.bash hook)"
            conda activate simulation
            rm -rf tskibd
            git clone git@github.com:bguo068/tskibd.git tskibd
            cd tskibd
            git checkout 8a3aba38067143bcc7934fb8d8a56124e7a88c92
            git submodule update --init --recursive
            meson build
            ninja -C build tskibd
            cd ../
            cp {curdir}/tskibd/build/tskibd {curdir}/bin/
            rm -rf tskibd
            """,
            shell=True,
            check=True,
            executable="bash",  # to sepcificy shell explicitly
        )
        print("install tskibd into bin/")
    else:
        print("tskibd already available in bin/ dir")


def install_hmmibd():
    curdir = Path(".").absolute()
    (curdir / "bin").mkdir(parents=True, exist_ok=True)
    if not (curdir / "bin/hmmIBD").exists():
        run(
            f"""
            eval "$(conda shell.bash hook)"
            conda activate simulation
            git clone https://github.com/glipsnort/hmmIBD.git
            cd hmmIBD
            git checkout a2f796ef8122d7f6b983ae9ac4c6fba35afcd3aa
            sed -i -e 's/const double rec_rate = 7.4e-7/const double rec_rate = 6.67e-7/' \
                    hmmIBD.c
            x86_64-conda_cos6-linux-gnu-gcc -o hmmIBD -O3 -lm -Wall hmmIBD.c
            cd ..
            cp hmmIBD/hmmIBD bin/hmmIBD
            rm -rf hmmIBD
             """,
            shell=True,
            check=True,
            executable="bash",
        )
        print("install hmmIBD into bin/")
    else:
        print("hmmIBD already available in bin/ dir")


def download_ibdne():
    if Path("bin/ibdne.jar").exists():
        print("ibdne already downloaded to bin/ibdne.jar")
    else:
        print("download ibdne.jar into bin/")
        url = "https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar"
        assert run(f"wget {url} -O bin/ibdne.jar", shell=True).returncode == 0


if __name__ == "__main__":

    if not is_conda_available():
        print("Conda is not installed. Please install Conda!", file=sys.stderr)
        sys.exit(-1)
    else:
        print("Conda is installed")

    if not is_simulation_env_installed():
        install_simulation_env("./env.yaml")
    else:
        print("'simulation' Conda environment is installed")

    install_tskibd()

    install_hmmibd()

    download_ibdne()
