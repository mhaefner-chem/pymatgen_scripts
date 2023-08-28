from ase.io import vasp
from gpaw import GPAW
from ase.parallel import parprint

structure = vasp.read_vasp('POSCAR')

magmoms = []
for i in structure.get_chemical_symbols():
    if "Fe" == i:
        magmoms.append(5.0)
    elif "V" == i:
        magmoms.append(2.0)
    else:
        magmoms.append(0.0)


structure.magmoms = magmoms

# Hubbard U
U = {'Co':':d,3.32','Cr':':d,3.7','Fe':':d,5.3','Mn':':d,3.9','Mo':':d,4.38','Ni':':d,6.2','V':':d,3.25','W':':d,6.2'}


calc = GPAW(mode='lcao',
            xc='PBE',
            basis='dzp',
            convergence = {'density':0.002},
            setups=U,
            kpts={'size':(2,2,2), 'gamma': True})
structure.calc = calc
energy = structure.get_potential_energy()
structure.calc.write('GPAW.gpw')
