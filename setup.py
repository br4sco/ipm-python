from setuptools import setup, find_packages

from ipm_python import __version__

setup(
    name="ipm_python",
    version=__version__,
    license="MIT License",
    keywords=["ipm", "furuta", "model"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "matplotlib",
        "scipy"
    ],
    packages=find_packages(
        where=".",
        include=["ipm_python"],
    ),
)
