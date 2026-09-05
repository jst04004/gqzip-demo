from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "_gqzip_cpp",
        [
            "src/bindings.cpp",
            "src/engine.cpp",
            "src/quantizer.cpp",
            "src/rans_codec.cpp",
            "src/simd_scanner.cpp",
        ],
        include_dirs=["include"],
        cxx_std=20,
    ),
]

setup(
    name="gqzip",
    version="1.0.0",
    description="Context-Adaptive Genomic Quality Score Quantization and Bounded-Memory Streaming FASTQ Compression Engine",
    author="Jonathan S. Taylor",
    author_email="contact@gqzip.org",
    url="https://github.com/jst04004/gqzip-demo",
    license="Apache-2.0",
    package_dir={"": "python"},
    packages=find_packages(where="python"),
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "psutil>=5.8.0",
        "pybind11>=2.10.0",
    ],
    entry_points={
        "console_scripts": [
            "gqzip=gqzip.cli:main",
        ],
    },
)
