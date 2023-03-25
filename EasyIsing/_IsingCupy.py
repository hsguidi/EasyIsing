import numpy as np
import cupy
from cupyx.scipy import ndimage as cuimage
from EasyIsing import Ising



        
class IsingCupy(Ising):
    def __init__(self, L: int, seed):
        super().__init__(L, seed)
        assert(L%2==0)
        self.name = 'Cupy'
        self.rng = cupy.random.default_rng(self.seed)
        self.lattice = self.rng.integers(2, size=(L,L), dtype=np.int8)*2-1
        self.mask = cupy.ones_like(self.lattice, dtype=np.int8)
        self.mask[::2, ::2] = 0
        self.mask[1::2, 1::2] = 0
        self.amask = (1-self.mask)
        self.kernel = cupy.array([[0,1,0],[1,0,1],[0,1,0]])
        
    def magnetization(self):
        return cupy.sum(self.lattice).get()

    def energy(self):
        k = cupy.array([[0,1,0],[1,0,0],[0,0,0]])
        n = cuimage.convolve(self.lattice, k, mode='wrap')
        return cupy.sum(-self.lattice * n).get()
    
    def update(self, mcs, temperature, field):
        beta=2.0/temperature
        for _ in range(mcs):         
            n = cuimage.convolve(self.lattice, self.kernel, mode='wrap')
            de = (field + n) *self.lattice
            self.lattice[self.rng.random(self.lattice.shape, dtype=np.float32) < cupy.exp(-beta * de, dtype=np.float32) * self.mask] *= -1

            n = cuimage.convolve(self.lattice, self.kernel, mode='wrap')
            de = (field + n) *self.lattice
            self.lattice[self.rng.random(self.lattice.shape, dtype=np.float32) < cupy.exp(-beta * de, dtype=np.float32) * self.amask] *= -1
        
        self.age += mcs
        self.E = self.energy()
        self.M = self.magnetization()
        

        
