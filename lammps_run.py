
# This script set up lammps run and retrieve atomization energies + coordinates
import os
import commands
import re

fileList = os.listdir('structures')
potEnergies = {} # potential energies
for file in fileList:
    strFile = os.path.join('structures',file)
    print 'Structure:', file
    fIn = open(strFile,'rU')
    linesIn =  fIn.readlines()
    fIn.close()

    noa =  len(linesIn) # number of atoms = number of lines (no blank lines)
    fOut = open('validation.str','w') #validation.str file
    linesOut =[]    

### validation.str    
    lineNum = 0
    count = 0 # count number of elements
    dict = {} # store element id (integer)
    for line in linesIn:
        lineNum += 1
        [element, pos_x, pos_y, pos_z] = line.split()
        if element not in dict.keys():
            count += 1
            dict[element] = count
        text = '%8d%8d%8d%16s%16s%16s' %(lineNum,dict[element],0,pos_x,pos_y,pos_z)
        linesOut.append(text)
    header = '''#%s

%d atoms

%d atom types
-50 50 xlo xhi
-50 50 ylo yhi
-50 50  zlo zhi

Atoms

''' %(file,noa,count)
    fOut.write(header)
    fOut.write('\n'.join(linesOut))
    fOut.close()

### in.reax
    fOut2 = open('in.reax','w')
    linesOut2 = []
    mass = { 'C':12.0107,'H':1.00794,'O':15.9994,'Al':26.982,'Cl':35.453 } # mass dictionary
    massLines = []
    elements = []
    for element in dict.keys():
        massLines.append('mass\t%d %f' %(dict[element],mass[element]))
        elements.append(element)
    text = '\n'.join(massLines)
    elements = ' '.join(elements)
    header2 = '''echo            both

units           real

atom_style      charge
read_data       validation.str
%s
pair_style      reax/c NULL

pair_coeff      * * ffield %s

compute reax all pair reax/c

variable eb      equal c_reax[1]
variable ea      equal c_reax[2]
variable elp     equal c_reax[3]
variable emol    equal c_reax[4]
variable ev      equal c_reax[5]
variable epen    equal c_reax[6]
variable ecoa    equal c_reax[7]
variable ehb     equal c_reax[8]
variable et      equal c_reax[9]
variable eco     equal c_reax[10]
variable ew      equal c_reax[11]
variable ep      equal c_reax[12]
variable efi     equal c_reax[13]
variable eqeq    equal c_reax[14]

neighbor        2.5 bin
neigh_modify    every 10 delay 0 check no

fix             2 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

thermo          1
thermo_style    custom step temp epair etotal press &
                v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa &
                v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq
# dump file 1
dump            1 all atom 1 dump.reaxc
# dump file 2
dump            4 all cfg 1 *.cfg mass type xs ys zs vx vy vz id q
dump_modify     4 element %s
# dump file 3
dump		5 all xyz 1 dump.xyz

minimize        1e-6 1e-8 10000 100000''' % (text,elements,elements)
    fOut2.write(header2)
    fOut2.close()

### lammps run
    cmd = './lmp_openmpi < in.reax'
    print "Command to run: ",cmd
    (status,output) = commands.getstatusoutput(cmd)

### read lammps log file
    fIn3 = open('log.lammps','rU')
    linesIn3 = fIn3.read()
    fIn3.close()
    match = re.search(r'\n\s*\d+\s+0\s+([-0-9.]+)\s+.*\nLoop\stime\sof\s',linesIn3)
    potEnergy = float(match.group(1))
    potEnergies[file] = potEnergy

print potEnergies

  
    
