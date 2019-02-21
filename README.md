# GeenuFF
Schema and API for a relational db that encodes gene models in an explicit, structured, and robust fashion.

## beta disclaimer

GeenuFF is currently _extremely_ beta and very unstable. 
We're keen to get feed back or ideas from the community
(even if it's just whether you think this could be useful to you
if developed further), but if you build on GeenuFF as it is now, 
you're doing so at your own risk.

## Motivation
We developed this to provide a way of unambiguously encoding gene models, 
(the way the DNA sequence is interpreted to produce proteins) that is both
robust to partial information and biological complexity.

A more extensive description can be found in `docs/full_description.md`

## Install
GeenuFF has been tested exclusively (so far) in python3.6

I would recommend installation in a virtual envinronment.
https://docs.python-guide.org/dev/virtualenvs/

The dustdas dependency is available from github.
From a directory of your choice (and preferably in a virtualenv):

```bash
git clone https://github.com/janinamass/dustdas.git
cd dustdas
python setup.py install
cd ..
```

Afterwards install GeenuFF in the same fashion.

```bash
git clone https://github.com/weberlab-hhu/GeenuFF.git
cd GeenuFF
pip install -r requirements.txt
python setup.py install
cd ..
```

And you might want to run the tests (sorry for the strict directory, will fix)
```bash
cd GeenuFF/geenuff
py.test
cd ../..
```

## Major plans
* Add a validation module to check structure of gene models.
* Add extraction of raw & mature transcript, CDS, and protein sequence as a demo application.
* Visualization.

## Thanks

To @janinamass for discussion and advice.