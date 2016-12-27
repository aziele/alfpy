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
------------

1. Python (https://www.python.org/) version 2.7 or >= 3.3
2. NumPy (http://www.numpy.org/).


Installation
------------

1. Recomended way to install is using `pip <https://pip.pypa.io/en/stable/installing/>`_
::

    pip install alfpy   # add --user if you don't have root


2. Alternatively, you can install from Github source code. Download and unzip it.
::

   cd alfpy
   python setup.py install


Alfpy usage
-----------

The examples of using Alfpy are available at: http://www.combio.pl/alfree/download/.


Testing
-------

To run tests, go to the alfpy source code directory and type::

    python -m unittest discover


If you want to test a specific file (e.g. ``test_word_distance.py``), type::

    python -m unittest tests.test_word_distance


License
-------

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
