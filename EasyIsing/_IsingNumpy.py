import numpy as np
from scipy import ndimage
from EasyIsing._Ising import Ising


class IsingNumpy(Ising):
    def __init__(self, L, seed):
        super().__init__(L, seed)
        self.name = 'Numpy'
        assert(L%2==0)
        self.mask = np.ones_like(self.lattice)
        self.mask[::2, ::2] = 0
        self.mask[1::2, 1::2] = 0
        self.amask = 1 - self.mask
        self.kernel = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])

    def update(self, mcs, temperature, field):
        beta=-2.0/temperature
        for _ in range(mcs):
            n = ndimage.convolve(self.lattice, self.kernel, mode='wrap')
            de = (field + n)*self.lattice
            self.lattice[self.rng.random(
                self.lattice.shape) < np.exp(beta * de) * self.mask] *= -1

            n = ndimage.convolve(self.lattice, self.kernel, mode='wrap')
            de = (field + n)*self.lattice
            self.lattice[self.rng.random(
                self.lattice.shape) < np.exp(beta * de) * self.amask] *= -1
        self.age += mcs
        self.E = self.energy()
        self.M = self.magnetization()

 
