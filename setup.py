from setuptools import setup, find_packages
from pip.req import parse_requirements

# parse_requirements() returns generator of pip.req.InstallRequirement objects
install_reqs = parse_requirements("requirements.txt", session=False)

# reqs is a list of requirement
reqs = [str(ir.req) for ir in install_reqs if ir.req is not None]

setup(name='igv_session',
      version='0.6.5',
      packages=find_packages(exclude=('tests*', 'docs', 'examples')),
      install_requires=reqs,
      entry_points={
          'console_scripts': [
              'igv_session = igv_session.igv_session:cli',
          ]
      }
      )
