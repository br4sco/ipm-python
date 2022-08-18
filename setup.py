from setuptools import setup, find_packages

__version__ = "0.0.4"

setup(
    name="ipm_furuta",
    version=__version__,
    license="MIT License",
    author="Oscar Eriksson",
    url="https://github.com/br4sco/ipm-python",
    keywords=["ipm", "furuta", "model"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["matplotlib", "numpy", "scipy>=1.4.0"],
    packages=find_packages(
        where=".",
        include=["ipm_furuta"],
    ),
)
