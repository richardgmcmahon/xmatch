from setuptools import setup, find_packages

setup(name='xmatch',
      version='1.0',
      license='GNU GPL v3+',
      description='object RA, Dec cross-matching using astropy',
      author='richardgmcmahon',
      author_email='rgm@ast.cam.ac.uk',
      packages=['xmatch',],
      install_requires=['astropy'],
      #packages=find_packages(include=['xmatch', 'xmatch.*'])
      #py_modules=['src','xmatch',]
     )
