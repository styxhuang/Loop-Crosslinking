# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 15:30:38 2018

@author: HuangMing

run Main.py -m DGEBA -mL 25 -c PACM -cL 15 -b 50 -Nm 4 -Nc 2 -r1 3 21 -r2 0 14 -cf 15. -B 1 -BT 2
"""
import argparse
import fileinput
import subprocess
import os
import sys
from shutil import copyfile
from shutil import move
import crosslinking
import procTOP
import readGRO
import readMol2

PACKMOL = 'packmol.inp'
SYSTEM = 'system'    
INI = 'init' 
BOND = 'bonded'
FF = 'oplsaa'
SOFT = ''
DATA = '$DATA'
GROMACS = 'gmx'
MPI = '' #mpirun -np 1

cwd = os.getcwd()
sys.path.append(cwd)

def Replace(filename, searchText, replaceText): #replace string in the file
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(searchText, replaceText), end='')
    
def MDSimulations(filename, cycle, monoNum, crosNum):
    command1 = 'obabel -i mol2 -o mol2 -fi {} -O {}.mol2 -h'.format(filename, INI)
    command2 = 'echo "" >> {}.mol2'.format(INI)
    command3 = '{}topolbuild -dir {} -ff {} -n {}'.format(SOFT, DATA, FF, INI)
    command4 = '{} editconf -f {}.gro -o box.gro -box 5 5 5'.format(GROMACS, INI)
    command5 = '{} grompp -f em.mdp -c box.gro -o min -maxwarn 10'.format(GROMACS)
    command6 = '{} {} mdrun -deffnm min -v'.format(MPI, GROMACS)
#    command7 = '{} grompp -f nvt.mdp -c min.gro -o nvt -maxwarn 10'.format(GROMACS)
#    command8 = '{} {} mdrun -deffnm nvt -v'.format(MPI, GROMACS)
    
    subprocess.call(command1, shell=True)
    subprocess.call(command2, shell=True)
    subprocess.call(command3, shell=True)
    
    fileName = '{}.top'.format(INI)
    Replace(fileName, 'ff{}'.format(INI), '{}.ff/forcefield'.format(FF))
    Replace(fileName, 'spce', '{}.ff/spc'.format(FF))
    Replace(fileName, 'ions', '{}.ff/ions'.format(FF))
    os.rename(fileName, 'topol.top')
    procTOP.Main('topol.top')
    os.mkdir('gmx')
    move('topol.top-end', 'gmx/topol.top')
    move('{}.gro'.format(INI), 'gmx/{}.gro'.format(INI))
    move('posre{}.itp'.format(INI), 'gmx/posre{}.itp'.format(INI))
    copyfile('../em.mdp', 'gmx/em.mdp')
    copyfile('../nvt.mdp', 'gmx/nvt.mdp')
    os.chdir('gmx')
    
    subprocess.call(command4, shell=True)
    subprocess.call(command5, shell=True)
    subprocess.call(command6, shell=True)
#    subprocess.call(command7, shell=True)
#    subprocess.call(command8, shell=True)
    
    os.mkdir('sim_result')
    move('min.gro', 'sim_result/min.gro')
    move('min.trr', 'sim_result/min.trr')
#    move('nvt.gro', 'sim_result/nvt.gro')
#    move('nvt.trr', 'sim_result/nvt.trr')
    move('topol.top', 'sim_result/topol.top')
    GRO2MOL2(cycle, monoNum, crosNum)

def GRO2MOL2(cycle, monoNum, crosNum):
    os.chdir('sim_result')
    top = 'topol.top'
    
    readGRO.Main('min', top, cycle)
#    readGRO.Main('nvt', top, cycle)
    command1 = 'obabel -i mol2 -o mol2 -fi min-stp-{}.mol2 -O min-dH.mol2 -d'.format(cycle)
    subprocess.call(command1, shell=True)
    
    newInfo = readMol2.InfoInput('min-dH.mol2', monLen, crosLen, monoNum, crosNum, dih=False) #TODO: After new molecule split method, sth happens here
#    readMol2.AtomInfoUpdate(oriInfo[0], newInfo[0], pos=True, charge=True)
    crosslinking.ExportMOL2(newInfo[5][0], 'min-dH.mol2'.format(cycle), newInfo[5][1:], newInfo[0], newInfo[1])

    copyfile('min-dH.mol2', '../../../tmp.mol2')
    os.chdir('../../')
#def GRO2MOL2(oriInfo, cycle, monLen, crosLen, update=True):
#    os.chdir('sim_result')
#    top = 'topol.top'
#    
#    readGRO.Main('min', top, cycle)
##    readGRO.Main('nvt', top, cycle)
#    command1 = 'obabel -i mol2 -o mol2 -fi min-stp-{}.mol2 -O min-dH.mol2 -d'.format(cycle)
#    subprocess.call(command1, shell=True)
#    
#    newInfo = readMol2.InfoInput('min-dH.mol2', monLen, crosLen, dih=False)
#    readMol2.AtomInfoUpdate(oriInfo[0], newInfo[0], pos=True, charge=True)
#    crosslinking.ExportMOL2(newInfo[5][0], 'min-dH.mol2'.format(cycle), newInfo[5][1:], oriInfo[0], oriInfo[1])
#
#    copyfile('min-dH.mol2', '../../../tmp.mol2')
#    os.chdir('../../')

def preProcess(fileName, cycle):
    MDSimulations(fileName, cycle)
    os.chdir('sim_result')
    top = 'topol.top'
    readGRO.Main('min', top, cycle)
    command1 = 'obabel -i mol2 -o mol2 -fi min-stp-{}.mol2 -O min-dH.mol2 -d'.format(cycle)
    subprocess.call(command1, shell=True)
    copyfile('min-dH.mol2', '../../../tmp.mol2')
    os.chdir('../../')
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Start Crosslink process')
    
    parser.add_argument('-m', '-mono', action = 'store', dest='mono',
                        help = 'input the monomer mol2 file', metavar = '')
    parser.add_argument('-mL', '-monoLen', action = 'store', dest='monoLen',
                        help = 'input the monomer atoms number w/o H', metavar = '')
    parser.add_argument('-c', '-cross', action = 'store', dest='cross',
                    help = 'input the crosslinker mol2 file', metavar = '')
    parser.add_argument('-cL', '-crossLen', action = 'store', dest='crossLen',
                    help = 'input the crosslinker atoms number w/o H', metavar = '')
    
    parser.add_argument('-b', '-box', dest='box',
                        default = 5, type = int,
                    help = 'input the box size(default box size is 5A)', metavar = '')
    parser.add_argument('-Nm', '-numberM', dest='num1',
                        default = 4, type = int,
                    help = 'input the monomer number(default monomer number is 4)', metavar = '')
    parser.add_argument('-Nc', '-numberC', dest='num2',
                        default = 2, type = int,
                    help = 'input the crosslinker number(default crosslinker number is 4)', metavar = '')
    parser.add_argument('-r1', '-reactM', action = 'append', dest='r1',
                        default = [], nargs = 2,
                    help = 'input the reactive atoms index in monomer, w/o hydrogen', metavar = '')
    parser.add_argument('-r2', '-reactC', action = 'append', dest='r2',
                        default = [], nargs = 2,
                    help = 'input the reactive atoms index in crosslinker w/o hydrogen', metavar = '')
    parser.add_argument('-cf', '-cutoff', dest='cutoff',
                        default = 5., type = float,
                    help = 'input the radius cutoff(default cutoff equals to 1nm)', metavar = '')
    parser.add_argument('-B', '-Bond', dest='bond',
                        default = 1, type = int,
                    help = 'input the bonds generation each cycle(default bond equals to 1)', metavar = '')
    parser.add_argument('-BT', '-bondTotal', dest='Bonds',
                        default = 1, type = int,
                    help = 'input the total bonds will be generated(default total bonds equals to 1', metavar = '')
       
    args = parser.parse_args()

#    parser.print_help() #For test only
    
    print('Monomer name \t\t\t= ', args.mono)
    print('Monomer atoms number(w/o H) \t= ', args.monoLen)
    print('Crosslinker name \t\t= ', args.cross)
    print('Crosslinker atoms number(w/o H)\t= ', args.crossLen)

    print('Box size \t\t\t= ', args.box)
    print('Monomer number \t\t\t= ', args.num1)
    print('Crosslinker number \t\t= ', args.num2)
    print('Reactive atoms on monomer \t= ', args.r1)
    print('Reactive atoms on crosslinker \t= ', args.r2)
    print('Cutoff \t\t\t\t= ', args.cutoff)
    print('Bonds generation each cycle \t= ', args.bond)
    print('Total bonds will be generated \t= ', args.Bonds)

##################################################################    
######Prepare packmol input file, update packmol.inp file
#    subprocess.Popen(command1, shell=True)
#    subprocess.Popen(command2, shell=True)
#    
    Replace(PACKMOL, 'FFF', SYSTEM) #TODO: remove comment later
    Replace(PACKMOL, 'MMM', args.mono)
    Replace(PACKMOL, 'CCC', args.cross)
    Replace(PACKMOL, 'SIZE', str(args.box))
    Replace(PACKMOL, 'NUM1', str(args.num1))
    Replace(PACKMOL, 'NUM2', str(args.num2))
    command1 = "obabel -i mol2 -o pdb -fi {}.mol2 -O {}.pdb".format(args.mono, args.mono)
    command2 = "obabel -i mol2 -o pdb -fi {}.mol2 -O {}.pdb".format(args.cross, args.cross)
    command3 = "packmol < packmol.inp > pack.log"
    command4 = "obabel -i pdb -o mol2 -fi {}.pdb -O {}.mol2 -d".format(SYSTEM, SYSTEM)
    command5 = "cp {}.mol2 tmp.mol2".format(SYSTEM)
    
    subprocess.call(command1, shell=True)
    subprocess.call(command2, shell=True)
    subprocess.call(command3, shell=True)
    subprocess.call(command4, shell=True)
    subprocess.call(command5, shell=True)
    
#    with open('run_1.sh', 'w') as out:
#        out.write(command1 + '\n')
#        out.write(command2 + '\n')
#        out.write(command3 + '\n')
#        out.write(command4 + '\n')
#        out.write(command5 + '\n')
#    subprocess.call('bash run_1.sh', shell=True) #call: process command until command finish, Popen, go on w/o waiting
#    subprocess.Popen('packmol < packmol.inp > pack.log', shell=True)
##################################################################

#    subprocess.Popen("obabel -i pdb -o mol2 -fi {}.pdb -O {}.mol2 -d".format(SYSTEM, SYSTEM), shell=True)

##### Start crosslinking loop #####
    bond_cycle = args.bond
    bonds_total = args.Bonds
    monoNum = args.num1
    crosNum = args.num2
    
    num = int(bonds_total/bond_cycle) + 1
    cycle = 0
    for i in range(0, num): #TODO: Problems happens, in the second loop, since bonds has been generated, the molecule length changed. Need to fix, doesnt change the molecule index when bonds generation
        if cycle <= num-1:
            cycle += 1
            print('cycle: ', cycle)
            os.mkdir('cycle{}'.format(cycle))
            os.chdir('cycle{}'.format(cycle))
            copyfile('../tmp.mol2','stp{}.mol2'.format(cycle))
            print('Cycle: ', cycle)
            fileName = 'stp{}.mol2'.format(cycle)
            outputName = '{}-{}.mol2'.format(BOND, cycle)
            monLen = int(args.monoLen)
            crosLen = int(args.crossLen)
            monR = [int(args.r1[0][0]), int(args.r1[0][1])]
            crosR = [int(args.r2[0][0]), int(args.r2[0][1])]
            cutoff = float(args.cutoff)
            bondsNum = int(args.bond)
            
#            if i == 0:
#                preProcess(fileName, cycle)
#            else:
            #Bond generation, info includes atomList, bondList
            oriInfo = crosslinking.main(fileName, outputName, monLen, crosLen, monR, crosR, cutoff, bondsNum, cycle, monoNum, crosNum)
        
            #mol2 --> gro file, excute md simulations
            MDSimulations(outputName, cycle, monoNum, crosNum)
        
            #convert gro --> mol2, info includes new atomsList, and bondsList
#            GRO2MOL2(oriInfo, cycle, monLen, crosLen)
            
            os.chdir('../')
        else:
            print("Loop has been finished, all possible bonds has been generated")
            sys.exit('All done!!')

    