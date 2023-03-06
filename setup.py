# needed for installation of refinegems
from setuptools import setup

setup(name='refineGEMs',
      version='1.0',
      description='Scripts to curate a GEM',
      author='Famke Baeuerle',
      author_email='famke.baeuerle@gmail.com',
      license='MIT',
      packages=['refinegems'],
      install_requires = [
            "cobra==0.22.0",
            "biopython==1.79",
            "bioservices",
            "memote==0.13.0",
            "pandas==1.2.4",
            "numpy==1.20.3",
            # "pyyaml==5.4.1",
            "gffutils==0.10.1",
            "markupsafe==2.0.1", 
            "depinfo==1.7.0",
            "sortedcontainers==2.4.0",
            "libchebipy==1.0.10",
            "ratelimit==2.2.1",
            "sqlalchemy==1.4.43",
            "venn==0.1.3",
            "ols-client==0.1.3",
            "seaborn==0.12.2",
            "sqlalchemy==1.4.43",
            "click==8.1.3"
            ],
      zip_safe=False)
