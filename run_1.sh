obabel -i mol2 -o pdb -fi None.mol2 -O None.pdb
obabel -i mol2 -o pdb -fi None.mol2 -O None.pdb
packmol < packmol.inp > pack.log
obabel -i pdb -o mol2 -fi system.pdb -O system.mol2 -d
cp system.mol2 tmp.mol2
