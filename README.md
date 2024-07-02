# genotype phenotype correlation case studies (gpc-cs)

This repository presents case studies using the [genophenocorr](https://github.com/monarch-initiative/genophenocorr) Python library
to explore genotype phenotype correlations in Mendelian disease.

## Set up

There are many ways of setting up Python projects and all of them should work here. One way of doing it is
with the virtual environment tool as follows. Note that we name our virtual environment gcpvenv, but you can name it
anything you want.

```
ython3 -m venv gpcvenv
source gpcvenv/bin/activate
pip install --upgrade pip
pip install genophenocorr ## do this in the genophenocorr directory until we can do PyPI
pip install jupyter ipykernel
python -m ipykernel install --name gpcvenv --user
jupyter-notebook
```
After this set the notebook kernel to gpcvenv, and the examples shown here should all work.
