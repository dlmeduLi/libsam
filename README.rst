Sam file related libraries
==========================

**Libsam** is a collection of `SAM/BAM
file <https://github.com/samtools/hts-specs>`__ related python library.
Packaged included in this library include:

-  samparser: A robust sam file format checker and parser

Install
=======

You can install this package using
`PyPi <https://pip.pypa.io/en/latest/index.html>`__:

.. code:: shell

    pip install libsam
    easy_install libsam

or you can **Download** the .tar.gz source file form
`PyPi <https://pypi.python.org/pypi>`__ and run the setup.py script:

.. code:: shell

    python setup.py install

Usage
=====

In your python script, add the following line(use *samparser* as an
example):

.. code:: python

    from libsam import samparser

then you can use the functions or classes in the corresponding
libraries:

.. code:: python

    sam = samparser.Sam()
    ...

Copyright
=========

Copyright (c) 2015 dlmeduLi@163.com
