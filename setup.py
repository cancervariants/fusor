"""Module for package and distribution"""
from setuptools import setup

exec(open("fusor/version.py").read())
setup(version=__version__)
