try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

def readme():
    with open('README.rst') as f:
        return f.read()

config = {
    'description': 'DEseq pipeline for differential DHS analysis',
    'author': 'Summer Elasady for GSK/Altius Institute',
    'download_url': 'https://github.com/Altius/deseq_pipeline/tarball/0.1',
    'author_email': 'summer@altius.org',
    'version': '0.1',
    'packages': ['deseq', 'normalization', 'tests'],
    'data_files': [('deseq', 'deseq.R')],
    'name': 'deseq_pipeline',
    'install_requires' : [
        'ConfigParser',
        'pandas',
        'argparse',
        'glob',
        'logging',
        'requests',
        'socket',
        'statsmodels'
      ],
    'include_package_data': 'True'
}

setup(**config)