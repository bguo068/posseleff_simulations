#! /usr/bin/env python3
import sys
from pathlib import Path
from subprocess import run


# Install tskibd
def install_tskibd():
    curdir = Path(".").absolute()
    (curdir / "bin").mkdir(parents=True, exist_ok=True)
    if not (curdir / "bin/tskibd").exists():
        run(
            f"""
            eval "$(pixi shell-hook)"
            rm -rf tskibd
            git clone https://github.com/bguo068/tskibd.git
            cd tskibd
            git checkout 8a3aba38067143bcc7934fb8d8a56124e7a88c92
            git submodule update --init --recursive
            # avoid error when compiled on macos
            mv tskit/c/VERSION tskit/c/VERSION.txt
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
            """
            eval "$(pixi shell-hook)"
            git clone https://github.com/glipsnort/hmmIBD.git
            cd hmmIBD
            git checkout a2f796ef8122d7f6b983ae9ac4c6fba35afcd3aa
            sed -i -e 's/const double rec_rate = 7.4e-7/const double rec_rate = 6.67e-7/' \
                    hmmIBD.c
            # avoid error when compiled on macos
            $CC -o hmmIBD -O3 -lm -Wall hmmIBD.c
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
        run(f"wget {url} --no-check-certificate -O bin/ibdne.jar", shell=True, check=True)


if __name__ == "__main__":

    install_tskibd()

    install_hmmibd()

    download_ibdne()
