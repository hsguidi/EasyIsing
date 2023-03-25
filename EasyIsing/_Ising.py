import numpy as np
from scipy import ndimage
from functools import cache


class Ising:
    """
    Define lattice,  basic methods and accumulate data.
    
    lattice : ndarray
        Contains the state of all sites.
        
    E : scalar
        Total interaction energy, $-\sum \sigma_i \sigma_j$.
        The value is update at end of update() method.
        
    M : scalar
        Total magnetization, $\sum \sigma_i$.
        The value is update at end of update() method.
        
    age : int
        Number of Monte Carlo steps updates since initialization.
        
    name : str
        Name of this instance.
        
    L : int
        Lattice side lenght
        
    L2 : int
        Lattice size, number of sites.
        
    seed : int
        Random number generator seed.
    
    rng : numpy.random.Generator
        Random number generator.
        
    """
    def __init__(self, L, seed):
        self.name = 'Python'
        self.age = 0
        self.L = L
        self.L2 = L * L
        self.seed = seed
        self.rng = np.random.default_rng(self.seed)
        self.lattice = self.rng.integers(2,
            size=self.L * self.L, dtype=np.int8).reshape(self.L, self.L) * 2 - 1
        self.E = 0
        self.M = 0

    def randomLattice(self):
        """
        Write random spin over lattice sites.
        """
        for i in range(self.L * self.L):
            self.lattice.flat[i] = self.rng.integers(2, dtype=np.int8) * 2 - 1
        self.E = self.energy()
        self.M = self.magnetization()

    def update(self, mcs, temperature, field):
        """
        Perform a Metropolis update on the lattice by executing 'mcs' Monte Carlo steps.
        The 'E' and 'M' variables are updated after 'mcs' steps.
        
        Parameters:
        -----------
        mcs : int
            The number of Monte Carlo Steps to be updated.
            Each Monte Carlo step updates L*L sites.
                
        temperature : scalar
            The reduced temperature.
                
        field : scalar
            The reduced magnetic field, given by H/J.
            
        Returns:
        --------
        None
        
        """
        @cache
        def SP(n,e):
            de = 2.0*n*e
            de += 2.0*field*e
            return np.exp(-de/temperature)
        
        for _ in range(mcs*self.L2):
            i,j = self.rng.integers(self.L, size=2)
            sp = SP(self.lattice[i,(j+1)%self.L] + self.lattice[(i+1)%self.L,j] + self.lattice[i,j-1] + self.lattice[i-1,j], self.lattice[i,j])
            if sp > self.rng.random():
                self.lattice[i,j] *= -1
                
        self.age += mcs
        self.E = self.energy()
        self.M = self.magnetization()

    def magnetization(self):
        """
        Calculate and return the total magnetization of the lattice,
        which is given by the sum of all the spin values.

        Returns
        -------
        M : float
            Total magnetization of the lattice.

        """
        return np.sum(self.lattice)

    def energy(self):
        """
        Calculate and return the total energy of the lattice, which is given by
        the sum of sigma_i * sigma_j over all pairs of interacting spins.

        Returns
        -------
        E : float
            Total interaction energy of the lattice.

        Notes
        -----
        The energy of the lattice represents the total amount of interaction
        between neighboring spins. In the Ising model, the energy is minimized
        when all spins are aligned, which corresponds to a ferromagnetic state.
        The energy is an important quantity in the Ising model, as it is used to
        calculate thermodynamic quantities such as the heat capacity and magnetic susceptibility.
        """

        k = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
        n = ndimage.convolve(self.lattice, k, mode='wrap')
        return -np.sum(self.lattice * n)

    def snapshot(self):
        """
        Map spin orientation -1 and 1 to 0 and 1.
        Return the lattice in compact bytearray format.
        
        Returns
        -------
        B : bytearray
            Lattice in bytearray format.

        
        """
        B = bytearray()
        for i in range(self.L2 // 8 + self.L2 % 8 != 0):
            s = ''.join([str((j + 1) // 2)
                        for j in self.lattice.flat[i:i + 8]])
            B += int(s, base=2).to_bytes(1, byteorder='big')
        return B

    def sampling(self, sampleSize, temperature, field, sampleStep=1):
        """
        Generate a sequence of lattice states by performing Metropolis
        updates and recording the energy and magnetization at each step.

        Parameters
        ----------
        sampleSize : int
            The number of samples to be generated.

        temperature : float
            The reduced temperature at which to simulate the Ising model.

        field : float
            The reduced magnetic field at which to simulate the Ising model.

        sampleStep : int, optional (default=1)
            The number of Monte Carlo Steps (MCS) updates between samples.
            As default, sampleStep=1 means that for every MCS, the energy
            and magnetization are calculated and collected for statistics.
                
        """
        E1, E2, M1, M2 = 0, 0, 0, 0
        for _ in range(sampleSize):
            self.update(sampleStep, temperature, field)
            E1 += self.E
            M1 += self.M
            E2 += self.E * self.E
            M2 += self.M * self.M
        N = sampleSize
        return {'temperature': temperature,
                'field': field,
                'age': self.age,
                'sampleSize': sampleSize,
                'length': self.L,
                'seed': self.seed,
                'energy1': E1 / N,
                'magnet1': M1 / N,
                'energy2': E2 / N,
                'magnet2': M2 / N}


