import io
import os
from os.path import dirname, join
from setuptools import setup, Command


# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def get_version(relpath):
  """Read version info from a file without importing it"""
  for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
    if "__version__" in line:
      if '"' in line:
        # __version__ = "0.9"
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]


setup(
    name='transcriptM',
    version=get_version("transcriptM/__init__.py"),
#    url='',
    license='GPL-3.0',
    author='Peter Sternes',
    author_email='peter.sternes@qut.edu.au',
    description='transcriptM - metatranscriptome pipeline with metagenome contamination correction',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['transcriptM'],
    package_data={'': [
            "transcriptM/*",
                       ]},
    data_files=[(".", ["README.md", "LICENSE.md"])],
    include_package_data=True,
    install_requires= [
    ],
    # install via conda: click, pandas, pyyaml, snakemake
    entry_points={
          'console_scripts': [
              'transcriptM = transcriptM.transcriptM:main'
          ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)

