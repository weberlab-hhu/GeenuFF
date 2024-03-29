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

A more extensive description can be found [here](https://weberlab-hhu.github.io/GeenuFF/).

## Install
GeenuFF has been tested in python3.8

I would recommend installation in a virtual environment.
https://docs.python-guide.org/dev/virtualenvs/

From a directory of your choice (and preferably in a virtualenv):

Clone and install GeenuFF:

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

Alternatively, install directly from github via pip:

```bash
pip install git+https://github.com/weberlab-hhu/GeenuFF.git
```
## usage
You can run `bash example.sh` for a quick start with public data.
 This will setup the folder 'three_algae', download public data in
 the expected format, and import it into a geenuff spec db for each
species. For more information please see 
[the api docs](https://weberlab-hhu.github.io/GeenuFF/api.html).

## Major plans
* Add a validation module to check structure of gene models.
* Add extraction of raw & mature transcript, CDS, and protein sequence as a demo application.
* Visualization.

## Thanks

To @janinamass for discussion and advice.
