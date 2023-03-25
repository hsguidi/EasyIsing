import numpy as np
from scipy import ndimage
from functools import cache
import EasyIsing


class IsingC(EasyIsing.Ising):
    """
    Define lattice, basic methods and load library 'EasyIsing.IsingC.so'.
    
    lattice : ndarray
        Contains the state of all sites.
        The address of this ndarray is passed to 'EasyIsing.IsingC.so'
        at initialization.
        You should not reallocate it.
    
    """
    def __init__(self, L: int, seed):
        super().__init__(L, seed)
        self.name = 'C'
        #!gcc -shared -Wl,-soname,ising -o IsingC.so -fPIC ising.c -lm
        self.machine = np.ctypeslib.load_library('EasyIsing.IsingC.so', '.')

        self.initMachine = self.machine.initSystem
        self.initMachine.restype = None
        self.initMachine.argtypes = [
            np.ctypeslib.ctypes.c_ssize_t,
            np.ctypeslib.ctypes.c_int,
            np.ctypeslib.ndpointer(dtype=np.int8, ndim=2, shape=(L, L), flags='C_CONTIGUOUS')]
        self.initMachine(self.L, self.seed, self.lattice)

        self.updateMachine = self.machine.update
        self.updateMachine.argtypes = [np.ctypeslib.ctypes.c_int,
                                       np.ctypeslib.ctypes.c_double,
                                       np.ctypeslib.ctypes.c_double]

    def update(self, mcs, temperature, field):
        deltaE = self.updateMachine(mcs, 1.0 / temperature, field)
        self.age += mcs
        self.E = self.energy()
        self.M = self.magnetization()


