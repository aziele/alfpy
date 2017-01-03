|Travis| |PyPI| |Landscape| |Codecov|

alfpy
=====

alfpy is a bionformatics Python package that provides alignment-free framework 
to compare biological sequences (DNA/RNA/protein) and infers their 
phylogenetic relationships. 

alfpy also contains Python scripts with user-friendly command-line interfaces 
that let you compare unaligned FASTA sequences with more than 40 distance methods.


Latest source code
------------------
The official source code repository is at: https://github.com/aziele/alfpy


Web sites
---------
alfpy is also available as a web app: http://www.combio.pl/alfree


Requirements
============

1. Python (https://www.python.org/) version 2.7 or >= 3.3
2. NumPy (http://www.numpy.org/).


Installation
============

Option 1: Get the latest official version
-----------------------------------------

Install the latest official version with `pip <https://pip.pypa.io/en/stable/installing/>`_
::

   sudo pip install alfpy

If you are not allowed to use `sudo`, install alfpy as user::

   sudo pip install --user alfpy



Option 2: Get the latest development version
--------------------------------------------

Get it using this shell command, which requires Git::

   git clone https://github.com/aziele/alfpy.git

If you don't feel like using git, just download the package manually as a `gzipped tarball <https://github.com/aziele/alfpy/archive/master.zip/>`_.

Unpack the zip package, go to the directory and run the installation::

   cd alfpy
   python setup.py install

or::

   python setup.py install --user

Alfpy usage
===========

The examples of using Alfpy are available at: http://www.combio.pl/alfree/download/.


Testing
=======

To run tests, go to the alfpy source code directory and type::

    python -m unittest discover


If you want to test a specific file (e.g. ``test_word_distance.py``), type::

    python -m unittest tests.test_word_distance


Contact
=======

Drop us any feedback at: bioinfo@amu.edu.pl or on twitter `@a_zielezinski <https://twitter.com/a_zielezinski>`_.

License
=======

alfpy is under the MIT license; see ``LICENSE.txt``. Distribution, 
modification and redistribution, incorporation into other software,
and pretty much everything else is allowed.


.. |Travis| image:: https://travis-ci.org/aziele/alfpy.svg?branch=master
    :target: https://travis-ci.org/aziele/alfpy


.. |PyPI| image:: https://img.shields.io/pypi/v/alfpy.svg?branch=master
    :target: https://pypi.python.org/pypi/alfpy

.. |Landscape| image:: https://landscape.io/github/aziele/alfpy/master/landscape.svg?style=flat
   :target: https://landscape.io/github/aziele/alfpy/master
   :alt: Code Health

.. |Codecov| image:: https://codecov.io/gh/aziele/alfpy/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/aziele/alfpy
