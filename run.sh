obabel -i mol2 -o pdb -fi DGEBA.mol2 -O DGEBA.pdb
obabel -i mol2 -o pdb -fi PACM.mol2 -O PACM.pdb
packmol < packmol.inp > pack.log
obabel -i pdb -o mol2 -fi system.pdb -O system.mol2 -d
