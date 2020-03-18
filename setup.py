from setuptools import setup

setup(
   name='geenuff',
   version='0.1',
   description='Schema and API for a relational db that encodes gene models in an explicit, structured, and robust fashion.',
   author='Alisandra Denton',
   packages=['geenuff', 'geenuff.base', 'geenuff.applications', 'geenuff.tests'],  #same as name
   install_requires=["intervaltree>=3.0.2", "SQLAlchemy>=1.3.12", "numpy>=1.18.1", "dustdas @ git+https://github.com/janinamass/dustdas@master", "pytest>=5.3.4"],
   zip_safe=False,
)
