from setuptools import setup


setup(
   name='alfpy',
   version='1.0.0',
   description="alignment-free package to compare DNA/RNA/protein sequences",
   author='Andrzej Zielezinski',
   keywords='alignment-free bioinformatics sequence DNA protein',
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
     'calc_word_bool2.py',
     'calc_word_comp.py',
     'calc_word_d2.py',
     'calc_word_ffp.py',
     'calc_word_rtd.py',
     'create_wordpattern.py'
   ],
   classifiers=[
     'License :: OSI Approved :: MIT License',
     'Programming Language :: Python :: 2',
     'Programming Language :: Python :: 2.7',
     'Programming Language :: Python :: 3',
     'Programming Language :: Python :: 3.3',
     'Programming Language :: Python :: 3.4',
     'Programming Language :: Python :: 3.5',
   ],   

)