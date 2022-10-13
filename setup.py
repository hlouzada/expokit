from pathlib import Path
from setuptools import find_packages
from numpy.distutils.extension import Extension
from numpy.distutils.core import setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="expokit",
    version="0.0.1",
    author="..",
    author_email="..",
    description="expokit package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux"
    ],
    ext_modules=[Extension(name = 'expokit.expokit',
                 sources = list(str(p) for p in Path('lib').glob('*.f')))],
    python_requires='>=3.6',
    zip_safe=False,
    install_requires=['numpy>=1.18']
)