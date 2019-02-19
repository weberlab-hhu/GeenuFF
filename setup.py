from setuptools import setup

setup(
   name='geenuff',
   version='0.1',
   description='Schema and API for a relational db that encodes gene models in an explicit, structured, and robust fashion.',
   author='Alisandra Denton',
   packages=['geenuff', 'geenuff.base', 'geenuff.applications', 'geenuff.tests'],  #same as name
)

