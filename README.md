
# ipm-python
A model of the Furuta pendulum based on the article [Furuta's Pendulum: A
Conservative Nonlinear Model for Theory and
Practise](https://doi.org/10.1155/2010/742894) by J. ́A. Acosta. The model
interface is based on [Interactive Programmatic
Modeling](https://dl.acm.org/doi/10.1145/3431387) by D. Broman.

Under most circumstances the imperative interface is preferred. The functional
style interface is there mainly for completeness.


Installing dependencies via the pip package manager:

```
pip install numpy scipy matplotlib
```

## Packaging

Make sure to run `make coherent` to ensure that the version from setup.py gets
propagated to the module.

Build the package using the setup.py script:

```
python setup.py sdist
```
