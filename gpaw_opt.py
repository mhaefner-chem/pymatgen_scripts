from ase.io import vasp
from gpaw import GPAW
from ase.optimize.bfgslinesearch import BFGSLineSearch

structure = vasp.read_vasp('POSCAR')
calc = GPAW(mode='lcao',
            xc='PBE',
            basis='dzp',
            convergence = {'density':0.002},
            kpts={'size':(4,4,4), 'gamma': True})
structure.calc = calc
energy = structure.get_potential_energy()

relax = BFGSLineSearch(atoms=structure, trajectory='bfgs_ls.traj', restart='bfgs_ls.pckl')
relax.run(fmax=0.05)
