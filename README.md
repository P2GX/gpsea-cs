# GPSEA case studies (gpsea-cs)

This repository presents case studies using the [GPSEA](https://github.com/p2gx/gpsea) Python library
to explore genotype phenotype correlations in Mendelian disease.

Each notebook presents one example of the application of GPSEA. The notebooks contain a section called **Summary** that
we used to generate the supplemental files for the GPSEA manuscript. There is not need to use the code in this section
for a new analysis. We strongly advise new users to follow the [GPSEA documentation](https://p2gx.github.io/gpsea/stable/index.html), 
which will be kept up to date with any further changes to the GPSEA code.

## Set up

There are many ways of setting up Python projects and all of them should work here. One way of doing it is
with the virtual environment tool as follows. Note that we name our virtual environment gcpvenv, but you can name it
anything you want.

```
python3 -m venv gpcvenv
source gpcvenv/bin/activate
pip install --upgrade pip
pip install gpsea
pip install jupyter ipykernel
python -m ipykernel install --name gpcvenv --user
jupyter-notebook
```
After this set the notebook kernel to gpcvenv, and the examples shown here should all work.

## gpseacs Python package

We created a Python package to streamline the creation of cohort summaries for the supplemental material.
Users of GPSEA will not need to use this package (unless they desire to create a similar supplemental file for
a new project).
