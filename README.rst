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

Use `pip <https://pip.pypa.io/en/stable/installing/>`_ to download, build and install alfpy and its dependencies::

    pip install alfpy


or download and unzip the source code, go to this directory at the command line, and type::

    sudo python setup.py install


Alfpy usage
-----------

The examples of using Alfpy are available at: http://www.combio.pl/alfree/download/.


Testing
-------

To run tests, go to the alfpy source code directory and type::

    python -m unittest discover


If you wanto to test a specific file (e.g. ``test_word_distance.py``), type::

    python -m unittest tests.test_word_distance


License
-------

alfpy is under the MIT license; see ``LICENSE.txt``. Distribution, 
modification and redistribution, incorporation into other software,
and pretty much everything else is allowed.