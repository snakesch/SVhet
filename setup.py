from setuptools import setup, find_packages
from setuptools.command.install import install
import compileall

class PostInstallCompile(install):
    """Post-installation script to compile Python files."""
    def run(self):
        install.run(self)
        compileall.compile_dir(self.install_lib, force=True)
        
setup(
    name="svhet",
    version="0.1.0",
    description="An accurate NGS-based structural variant filtering tool using heterozygous sites",
    packages=find_packages(where="core"),
    package_dir={"": "core"},
    cmdclass={'install': PostInstallCompile},
    install_requires=[
        "cyvcf2>=0.31.0",
        "pysam>=0.23.0",
        "pybedtools>=0.12.0",
        "tqdm>=4.67.1",
        "numpy>=1.26.4,<2",
        "scipy>=1.15.2"
    ],
    entry_points={
        "console_scripts": [
            "svhet=svhet.cli:main",
        ],
    },
    author="Louis She",
    author_email="snakesch@connect.hku.hk",
    keywords="genomics, structural variants, variant detection",
    url="https://github.com/snakesch/svhet",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bioinformatics",
    ],
    python_requires=">=3.11",
)
