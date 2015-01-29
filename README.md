# Sam file related libraries

**Libsam** is a collection of [SAM/BAM file](https://github.com/samtools/hts-specs) related python library.
Packaged included in this library include:

 - samparser: A robust sam file format checker and parser

Install
===
You can install this package using [PyPi](https://pip.pypa.io/en/latest/index.html):
```shell
pip install libsam
easy_install libsam
```
or you can **Download** the .tar.gz source file form [PyPi](https://pypi.python.org/pypi) and run the setup.py script:
```shell
python setup.py install
```
Usage
==
In your python script, add the following line(use *samparser* as an example):
```python
from libsam import samparser
```
then you can use the functions or classes in the corresponding libraries:
```python
sam = samparser.Sam()
...
```

Copyright
==
Copyright (c) 2015 [dlmeduLi@163.com](mailto:dlmeduLi@163.com)
