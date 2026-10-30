#! /usr/bin/env python3
import shutil
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
            if [ "$(uname)" = "Darwin" ]; then export LDFLAGS="-lSystem"; fi
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
            if [ "$(uname)" = "Darwin" ]; then export LDFLAGS="-lSystem"; fi
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


# Build SLiM where conda has no build of it (linux-aarch64)
def install_slim():
    curdir = Path(".").absolute()
    (curdir / "bin").mkdir(parents=True, exist_ok=True)
    if shutil.which("slim") or (curdir / "bin/slim").exists():
        print("slim already available")
        return
    run(
        f"""
        set -e
        eval "$(pixi shell-hook)"
        rm -rf SLiM
        git clone --depth 1 --branch v4.1 https://github.com/MesserLab/SLiM.git
        cmake -G Ninja -S SLiM -B SLiM/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
        cmake --build SLiM/build --target slim -j 4
        cp SLiM/build/slim {curdir}/bin/slim
        rm -rf SLiM
        """,
        shell=True,
        check=True,
        executable="bash",
    )
    print("install slim into bin/")


if __name__ == "__main__":

    install_slim()

    install_tskibd()

    install_hmmibd()

    download_ibdne()
