# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.md') as f:
        readme = f.read()

with open('LICENSE') as f:
    license = f.read()
    
setup(
    name='CAIROpy',
    version='0.1.1',
    description='python interface package for CAIRO software',
    long_description=readme,
    classifiers=[
        'Development status :: 1 - Alpha',
        'License :: CC-By-SA2.0',
        'Programming Language :: Python',
        'Topic :: Data Analysis'
    ],
    author='Antoine Marchal',
    author_email='antoine.marchal@anu.edu.au',
    url='https://github.com/antoinemarchal/CAIRO',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=[
            'numpy',
            'matplotlib'
    ],
    include_package_data=True
)
