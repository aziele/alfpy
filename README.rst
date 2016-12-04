alfpy
=====

alfpy is a bionformatics Python package that provides alignment-free framework 
to compare biological sequences (DNA/RNA/protein) and infers their 
phylogenetic relationships. 

alfpy also contains Python scripts with user-friendly command-line interfaces 
that let you compare unaligned FASTA sequences with more than 40 distance methods.

The official source code repository is at: https://github.com/aziele/alfpy

alfpy is also available as a web app: http://www.combio.pl/alfree


Requirements
------------

alfpy requires:

1. Python (https://www.python.org/) version 2.7 or >= 3.3
2. NumPy (http://www.numpy.org/).


Installation
============

Install using pip
-----------------

Use `pip<http://www.pip-installer.org/>`_ to download, build, and install alfpy and its dependencies:
    pip install alfpy


Install from GitHub
-------------------

You may have to build and install alfpy. Download and unzip the
source code, go to this directory at the command line, and type::

    sudo python setup.py install

Here you can replace ``python`` with a specific version, e.g. ``python3.5``.


Alfpy usage
===========

The examples of using Alfpy are available at: http://www.combio.pl/alfree/download/.


License
-------

alfpy is under the MIT license; see doc/LICENSE.txt. Distribution, 
modification and redistribution, incorporation into other software, and 
pretty much everything else is allowed.