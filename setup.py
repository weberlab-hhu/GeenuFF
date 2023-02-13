try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
   name='geenuff',
   version='0.3.1',
   description='Schema and API for a relational db that encodes gene models in an explicit, structured, and robust fashion.',
   author='Alisandra Denton, Felix Stiehler',
   packages=['geenuff', 'geenuff.base', 'geenuff.applications', 'geenuff.applications.exporters', 'geenuff.tests'],  #same as name
   package_data={'geenuff': ['testdata/*.fa', 'testdata/*.gff3', 'testdata/*.gff']},
   install_requires=["intervaltree>=3.0.2", "SQLAlchemy>=1.3.12", "numpy>=1.18.1", "dustdas @ git+https://github.com/janinamass/dustdas@master", "pytest>=5.3.4"],
   scripts=['import2geenuff.py'],
   zip_safe=False,
)
