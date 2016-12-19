from setuptools import setup

# Read a __version__
exec(open('alfpy/version.py').read())

setup(
   name='alfpy',
   version=__version__,
   description="Alignment-free package to compare DNA/RNA/protein sequences (bioinformatics).",
   author='Andrzej Zielezinski',
   keywords='alignment-free bioinformatics sequence DNA protein homology phylogeny',
   license="MIT",
   author_email='andrzejz@amu.edu.pl',
   url="http://www.combio.pl/alfree",
   packages=['alfpy', 'alfpy.utils', 'alfpy.utils.data'],
   #setup_requires=["numpy"],
   install_requires=["numpy"],
   scripts=[
     'calc_bbc.py',
     'calc_dna2d.py',
     'calc_fcgr.py',
     'calc_lempelziv.py',
     'calc_ncd.py',
     'calc_wmetric.py',
     'calc_word.py',
     'calc_word_bool.py',
     'calc_word_sets.py',
     'calc_word_comp.py',
     'calc_word_d2.py',
     'calc_word_ffp.py',
     'calc_word_rtd.py',
     'create_wordpattern.py'
   ],
   classifiers=[
     'License :: OSI Approved :: MIT License',
     'Environment :: Console',
     'Operating System :: MacOS',
     'Operating System :: Microsoft :: Windows',
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