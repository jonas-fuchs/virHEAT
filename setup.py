from setuptools import setup, find_packages
from virheat import __version__, _program

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='virheat',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=__version__,
    python_requires=">=3.9",
    license_files=('LICENSE'),
    packages=find_packages(),
    install_requires=[
        "matplotlib>=3.5.1",
        "numpy>=1.23.3",
    ],
    description='virHEAT creates a heatmap from vcf files and maps positions onto a reference genome.',
    url='https://github.com/jonas-fuchs/virHEAT',
    author='Dr. Jonas Fuchs',
    author_email='jonas.fuchs@uniklinik-freiburg.de',
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)"
    ],
    entry_points="""
    [console_scripts]
    {program} = virheat.command:main
    """.format(program=_program),
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
