"""
EasyIsing
=====

EasyIsing is a Python package that provides a simple Ising model simulation on a squared lattice with the Metropolis algorithm.

Available subpackages
---------------------

EasyIsing.Ising:
    This is the base class that provides a standard Metropolis implementation of the Ising model on a squared lattice.
    
EasyIsing.IsingC:
    This is a subclass of EasyIsing.Ising where the EasyIsing.IsingC.update() method calls the static library EasyIsing.IsingC.so to perform the Metropolis update.
    
EasyIsing.IsingNumpy:
    This is a subclass of EasyIsing.Ising where the spin selection in the Metropolis algorithm is parallelized by a checkerboard scheme.

EasyIsing.IsingCupy:
    This is the same as EasyIsing.IsingNumpy, but the calculation is done by a GPU.


Viewing documentation using IPython
-----------------------------------

Start IPython and import `EasyIsing`.
To see which functions are available in `EasyIsing`,
type ``EasyIsing.Ising.<TAB>`` (where ``<TAB>`` refers to the TAB key).
To view the docstring for a function, use
``EasyIsing.Ising?<ENTER>`` to view the docstring and ``EasyIsing.Ising??<ENTER>``
to view the source code.

Lattice in-place operation
-----------------------------
The variable 'lattice' should not be reallocated, especially in the 'EasyIsing.IsingC' subclass, where the 'lattice' address is passed to 'EasyIsing.IsingC.so' only at initialization.

"""

from EasyIsing._Ising import Ising
from EasyIsing._IsingNumpy import IsingNumpy


import importlib.util
if importlib.util.find_spec("cupy") is not None:
	from EasyIsing._IsingCupy import IsingCupy
else:
    print('Cupy not found.')
	
	
from EasyIsing._IsingC import IsingC
from os.path import exists
from os import system

#!gcc -shared -Wl,-soname,ising -o IsingC.so -fPIC ising.c -lm
sopath = 'EasyIsing.IsingC.so'
if not exists(sopath):
	srcpath = __path__[0] + '/ising.c'
	system('gcc -shared -Wl,-soname,ising -o '+sopath+' -fPIC '+srcpath+' -lm')
