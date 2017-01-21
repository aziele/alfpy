from setuptools import setup

# Read a __version__
exec(open('alfpy/version.py').read())

# Long description
fh = open('README.rst', encoding='utf-8')
long_description = fh.read()
fh.close()

setup(
   name='alfpy',
   version=__version__,
   description="Alignment-free package to compare DNA/RNA/protein sequences (bioinformatics).",
   long_description=long_description,
   author='Andrzej Zielezinski',
   keywords='alignment-free bioinformatics sequence DNA protein homology phylogeny',
   license="MIT",
   author_email='andrzejz@amu.edu.pl',
   url="http://www.combio.pl/alfree",
   packages=['alfpy', 'alfpy.utils', 'alfpy.utils.data'],
   #setup_requires=["numpy"],
   install_requires=["numpy"],
   scripts=[
     'bin/calc_bbc.py',
     'bin/calc_graphdna.py',
     'bin/calc_fcgr.py',
     'bin/calc_lempelziv.py',
     'bin/calc_ncd.py',
     'bin/calc_wmetric.py',
     'bin/calc_word.py',
     'bin/calc_word_bool.py',
     'bin/calc_word_sets.py',
     'bin/calc_word_cv.py',
     'bin/calc_word_d2.py',
     'bin/calc_word_ffp.py',
     'bin/calc_word_rtd.py',
     'bin/create_wordpattern.py'
   ],
   classifiers=[
     'License :: OSI Approved :: MIT License',
     'Environment :: Console',
     'Operating System :: MacOS',
     'Operating System :: POSIX :: Linux',     
     'Programming Language :: Python :: 2',
     'Programming Language :: Python :: 2.7',
     'Programming Language :: Python :: 3',
     'Programming Language :: Python :: 3.3',
     'Programming Language :: Python :: 3.4',
     'Programming Language :: Python :: 3.5',
     'Topic :: Scientific/Engineering',
     'Topic :: Scientific/Engineering :: Bio-Informatics',     
   ],   

)