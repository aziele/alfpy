from alfpy import word_vector

fh = open('README.rst')


try:
    word_vector.read_weightfile(fh)
except Exception:
    pass
fh.close()