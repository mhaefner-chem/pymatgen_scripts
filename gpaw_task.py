from ase.io import vasp
from gpaw import GPAW
from ase.optimize import LBFGS

structure = vasp.read_vasp('POSCAR')
calc = GPAW(mode='lcao',
            xc='PBE',
            txt='GPAW.out',
            basis='dzp',
            kpts={'size':(4,4,4), 'gamma': True})
structure.calc = calc
energy = structure.get_potential_energy()

relax = LBFGS(atoms=structure, trajectory='lbfgs.traj', restart='lbfgs.pckl')
relax.run(fmax=0.05)

#calc.write('GPAW.gpw')
print(f' Energy: {energy:5.2f} eV')
