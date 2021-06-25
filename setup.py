from setuptools import setup, find_packages
import numpy

setup(
    name='papakuhi',
    version='0.2',
    include_package_data = True,
    include_dirs=[numpy.get_include()],
    packages=find_packages()
)
