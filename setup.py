# needed for installation of refinegems
from setuptools import setup

setup(name='refinegems',
      version='1.2',
      description='Scripts to curate a GEM',
      author='Famke Baeuerle',
      author_email='famke.baeuerle@gmail.com',
      license='MIT',
      packages=['refinegems'],
      install_requires = ["cobra==0.22.0",
            "biopython==1.79",
            "bioservices==1.7.11",
            "memote==0.11.1",
            "escher==1.7.3",
            "pandas==1.2.4",
            "numpy==1.20.3",
            "pyyaml==5.4.1",
            "psycopg2-binary==2.9.1",
            "gffutils==0.10.1"],
      zip_safe=False)

# switching to setuptools 57.5.0 fixed the use_2to3 error of installation of suds_jurko
# help(function) shows docstring of said function (use q to quit)
# python3 -m build
# twine upload --repository testpypi dist/*
# always change version before uploading / remove old dist