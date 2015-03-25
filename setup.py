import os
from setuptools import setup, Command

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "libsam",
    version = "0.1.8",
    author = "dlmeduLi",
    author_email = "dlmeduLi@163.com",
    description = ("Bio-Informatics sam file libraries."),
    license = "BSD",
    keywords = "bioinformatic samfile parser",
    url = "https://github.com/dlmeduLi/libsam",
    packages=['libsam'],
    long_description=read('README.rst'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
    ],
    # ... Other setup options
    cmdclass={
        'clean': CleanCommand,
    }
)