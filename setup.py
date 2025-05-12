from setuptools import setup, find_packages

setup(
    name='grale',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'astropy',
        'matplotlib'
    ],
)
