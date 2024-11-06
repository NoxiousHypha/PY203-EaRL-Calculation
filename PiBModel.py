import numpy as np
import scipy.constants as con

# This is a simple script to run a 2D particle-in-a-box simulation of
# multi-ringed hydrocarbons (napthelene, anthracene, tetracene)
# for comparison against lab-reported absorption spectra

class HCRing:                       # We
    def __init__(self, n):
        self.dbonds = 2 * n + 1     # Double-bonds of linear rings
        self.homo = self.dbonds     # Highest occupied molecular orbital
        self.lumo = self.homo + 1   # Lowest unoccupied molecular orbital
        self.ly = 0.28e-9           # The height of a chain of aromatic rings
        self.lx = (0.24e-9) * n     # Chain width


def listEnergies(molecule):
    E0 = (con.hbar**2 * con.pi**2)/(2*con.m_e)
    earray = []
    nxarray = []
    nyarray = []
    for nx in range(1,molecule.lumo + 3):
        for ny in range(1,4):
            energy = (E0 * ((nx**2 / (molecule.lx)**2) + (ny**2 / (molecule.ly)**2)))/con.e
            earray.append([energy])
            nxarray.append([nx])
            nyarray.append([ny])
    a = np.hstack([nxarray,nyarray,earray])
    sort = a[a[:,2].argsort()]
    return sort

naphthelene = HCRing(2)
anthracene = HCRing(3)
tetracene = HCRing(4)

nEnergies = listEnergies(naphthelene)
aEnergies = listEnergies(anthracene)
tEnergies = listEnergies(tetracene)

np.savetxt('./naphtheleneEnergies.txt',nEnergies, delimiter=', ', header='Napthelene Calculated Energies, (nx, ny, energy)')
np.savetxt('anthraceneEnergies.txt',aEnergies, delimiter=', ', header='Anthracene Calculated Energies, (nx, ny, energy)')
np.savetxt('tetraceneEnergies.txt',tEnergies, delimiter=', ', header='Tetracene Calculated Energies, (nx, ny, energy)')