
# This script set up lammps run and retrieve atomization energies + coordinates
import os
import commands
import re
import sys
#import subprocess

######

def clean():
    cmdList = ['rm *.str','rm in.*','rm -f log.lammps log.cite']
    for cmd in cmdList:
        (status,output) = commands.getstatusoutput(cmd)
#        subprocess.check_call(cmd)
    
######    

def generate_inputs():
# generate inputs for lammps run
    fileIn1 = open('geo','rU')
    inLines1 = fileIn1.readlines()
    fileIn1.close()

    mass = { 'C':12.0107,'H':1.00794,'O':15.9994,'Al':26.982,'Cl':35.453 } # mass dictionary

    structures = []
    for line in inLines1:
        keyword = line.split()
        if keyword:
            if keyword[0] == 'DESCRP':
                cName = keyword[1] # compond/species name
                structures.append(cName)
                #print cName
                elementDict = {}
                lines = []
                elementList = []
                massList = []
                filename = cName + '.str'
                filename2= 'in.' + cName
                strFile = open(filename,'w') # open str file for writing
                noe = 0 #number of element types
                noa = 0
            if keyword[0] == 'HETATM':
                noa += 1
                element = keyword[2]
                if element not in elementDict.keys():
                    noe += 1
                    elementDict[element] = noe
                    elementList.append(element)
                    massList.append('mass\t%d %f' %(elementDict[element],mass[element]))
                eid = elementDict[element] #element id
                lines.append('%8s%8d%8d%16s%16s%16s' % (keyword[1],eid,0,keyword[3],keyword[4],keyword[5]))
            if keyword[0] == 'END':
                header = '''#%s

%d atoms

%d atom types
-50 50 xlo xhi
-50 50 ylo yhi
-50 50  zlo zhi

Atoms

''' %(cName,noa,noe)
                strFile.write(header)
                strFile.write('\n'.join(lines))
                strFile.close() # 'validation.str' file generated
                
                inFile = open(filename2,'w') # input script
                elementList = ' '.join(elementList)
                massList = '\n'.join(massList)
                if cName == 'C3H8_ob':
                    fixstyle = '''group carbon type 1
fix 99 carbon setforce 0.0 0.0 0.0
'''
                else:
                    fixstyle=''

                content = '''echo            both

units           real

atom_style      charge
read_data       %s
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
%s

thermo          1
thermo_style    custom step temp epair etotal press &
                v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa &
                v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq
# dump file 1
#dump            1 all atom 1 dump.reaxc
# dump file 2
# dump file 3
dump		5 all xyz 1 dump.xyz
# dump file 4
dump        10 all custom 1 dump.custom id element q x y z
dump_modify     10 element %s

minimize        1e-6 1e-8 10000 100000''' % (filename,massList,elementList,fixstyle,elementList)
                inFile.write(content)
                inFile.close()

    return structures

######

def trainset_in():
# read information in trainset.in

    f = open('trainset.in')
    lines = f.readlines()
    f.close()
    chargeList = []
    geometryList = []
    energyList = []
    #CHARGE
    for line in lines:
        keywords = line.split()
        if keywords == ['CHARGE']:
            status = 'CHARGE'
        if keywords == ['GEOMETRY']: status = 'GEOMETRY'
        if keywords == ['ENERGY']: status = 'ENERGY'
        if ((keywords) and (keywords[0][0] != '#') and (line[:3] not in ['END','CHA','GEO','ENE'])):
            if status == 'CHARGE':
                struc = keywords[0]
                weight = float(keywords[1])
                atomID = keywords[2]
                litVal = float(keywords[3])
                chargeList.append({'struc':struc,'weight':weight,'atomID':atomID,'litVal':litVal})
            if status == 'GEOMETRY':
                struc = keywords[0]
                weight = float(keywords[1])
                litVal = float(keywords[-1])
                atoms = keywords[2:-1]
                geometryList.append({'struc':struc,'weight':weight,'atoms':atoms,'litVal':litVal})
            if status == 'ENERGY':
                weight = float(keywords[0])
                [struc, divider] = keywords[2].split('/')
                divider = float(divider)
                reaction = keywords[1:-1]
                litVal = float(keywords[-1])
                energyList.append({'struc':struc,'weight':weight,'litVal':litVal,'reaction':reaction,'divider':divider})
    trainset_inDict = {'CHARGE':chargeList,'GEOMETRY':geometryList,'ENERGY':energyList}
    return trainset_inDict

 
 

######

def reax_run(list):
    reaxDict = {}
    pEnergyDict = {}
    for struc in list:
        strucDict= {}
#        baseEnergy = 0
        #print struc
        cmd = './lmp_openmpi < in.%s' %struc
        (status,output) = commands.getstatusoutput(cmd)  #lammps run
        #print 'lammps run status',struc,':',status
        if status!=0: print output
#        subprocess.check_call(cmd)

        f =  open('log.lammps','r')
        content = f.read()
        f.close()
        match = re.search(r'\n\s*\d+\s+0\s+([-0-9.]+)\s+.*\nLoop\stime\sof\s',content)
        pEnergyDict[struc] = float(match.group(1)) #potential energy

        #read dump.custom file
        f = open('dump.custom','r')
        content = f.readlines()
        f.close()
        noa = int(content[3]) # number of atoms
        for line in content[-noa:]:
            keywords = line.split()
            atomDict = {'element':keywords[1],'q':float(keywords[2]),'coordinates':[float(keywords[3]),float(keywords[4]),float(keywords[5])]}
            strucDict[keywords[0]] = atomDict
#            baseEnergy += pEnergyDict[atomDict['element']]
#        aEnergy = pEnergyDict[struc] - baseEnergy # Atomization Energy
#        strucDict['aEnergy'] = aEnergy
        reaxDict[struc] = strucDict

    for struc in reaxDict.keys(): #calculate atomization energy
        baseEnergy = 0
        strucDict = reaxDict[struc]
        for atomID in strucDict.keys():
            baseEnergy += pEnergyDict[strucDict[atomID]['element']]
        aEnergy = pEnergyDict[struc] - baseEnergy
        reaxDict[struc]['aEnergy'] = aEnergy

###########
#    print '===reaxDict check'
#    for struc in reaxDict.keys():
#        print struc
#        strucDict = reaxDict[struc]
#        for key in strucDict.keys():
#            print key,strucDict[key]
#    print '=====reaxDict check end'

    return reaxDict
            


######

def sortFun(list):
    return list[-1]

######
def calculate_error(reaxDict,trainset_inDict):
    import math

    f = open('log.err','w')
    errDict = {}
    totalErr = 0
    chargeErr = []

    f.write('CHARGE\n')
    f.write('structure\tatomID\treaxVal\tlitVal\tweight\terror\ttotalError\n')
    for dict in trainset_inDict['CHARGE']:
        struc = dict['struc']
        weight = dict['weight']
        atomID = dict['atomID']
        litVal = dict['litVal']
        reaxVal = reaxDict[struc][atomID]['q']
        elements = reaxDict[struc][atomID]['element']
        err = ((reaxVal-litVal)/weight)**2.0
        chargeErr.append([struc,weight,atomID,reaxVal,litVal,elements,err])
#    totalErr = 0
#    chargeErr = sorted(chargeErr,key = sortFun) #sorting
    for i in range(len(chargeErr)):
        totalErr += chargeErr[i][-1]
        chargeErr[i].append(totalErr)
        line = chargeErr[i]
        f.write('%s\t%s/%s\t%f\t%f\t%f\t%f\t%f\n' %(line[0],line[2],line[5],line[3],line[4],line[1],line[6],line[7]))

    bondErr = []
    angleErr = []
    for dict in trainset_inDict['GEOMETRY']:
        struc = dict['struc']
        weight = dict['weight']
        atoms = dict['atoms']
        litVal = dict['litVal']

        nop = len(atoms) # number of atoms
        if nop == 2 : # bond
            atomID = atoms[0]
            coord1 = reaxDict[struc][atomID]['coordinates']
            atomID = atoms[1]
#####
#            print reaxDict[struc].keys()
#            print struc,atomID
            coord2 = reaxDict[struc][atomID]['coordinates']
            bondLen = ((coord1[0]-coord2[0])**2.0+(coord1[1]-coord2[1])**2.0+(coord1[2]-coord2[2])**2.0)**0.5
            reaxVal = bondLen
            elements = [reaxDict[struc][atoms[0]]['element'],reaxDict[struc][atoms[1]]['element']]
            err = ((reaxVal-litVal)/weight)**2.0
            bondErr.append([struc,weight,atoms,reaxVal,litVal,elements,err])

        if nop == 3 : # angle
            atomID = atoms[0]
            coord1 = reaxDict[struc][atomID]['coordinates']
            atomID = atoms[1]
            coord2 = reaxDict[struc][atomID]['coordinates']
            atomID = atoms[2]
            coord3 = reaxDict[struc][atomID]['coordinates']

            bondVec1 = [coord1[0]-coord2[0],coord1[1]-coord2[1],coord1[2]-coord2[2]]
            bondVec2 = [coord3[0]-coord2[0],coord3[1]-coord2[1],coord3[2]-coord2[2]]
            bondLen1 = ((bondVec1[0])**2.0+(bondVec1[1])**2.0+(bondVec1[2])**2.0)**0.5
            bondLen2 = ((bondVec2[0])**2.0+(bondVec2[1])**2.0+(bondVec2[2])**2.0)**0.5

            dotProd = bondVec1[0]*bondVec2[0]+bondVec1[1]*bondVec2[1]+bondVec1[2]*bondVec2[2]
            valAngle = math.degrees(math.acos(dotProd/(bondLen1*bondLen2)))
            elements = [reaxDict[struc][atoms[0]]['element'],reaxDict[struc][atoms[1]]['element'],reaxDict[struc][atoms[2]]['element']]
            reaxVal = valAngle
            err = ((reaxVal-litVal)/weight)**2.0
            angleErr.append([struc,weight,atoms,reaxVal,litVal,elements,err])
    
    f.write('BONDLENGTH\n')
    f.write('structure\tatom1\tatom2\treaxVal\tlitVal\tweight\terror\ttotalError\n')
#    bondErr = sorted(bondErr,key = sortFun) #sorting
    for i in range(len(bondErr)):
        totalErr += bondErr[i][-1]
        bondErr[i].append(totalErr)
        line = bondErr[i]
        f.write('%s\t%s/%s\t%s/%s\t%f\t%f\t%f\t%f\t%f\n' %(line[0],line[2][0],line[5][0],line[2][1],line[5][1],line[3],line[4],line[1],line[6],line[7]))

    f.write('VALANCEANGLE\n')
    f.write('structure\tatom1\tatom2\tatom3\treaxVal\tlitVal\tweight\terror\ttotalError\n')
#    angleErr = sorted(angleErr,key = sortFun) #sorting
    for i in range(len(angleErr)):
        totalErr += angleErr[i][-1]
        angleErr[i].append(totalErr)
        line = angleErr[i]
        f.write('%s\t%s/%s\t%s/%s\t%s/%s\t%f\t%f\t%f\t%f\t%f\n' %(line[0],line[2][0],line[5][0],line[2][1],line[5][1],line[2][2],line[5][2],line[3],line[4],line[1],line[6],line[7]))

    energyErr = []
    for dict in trainset_inDict['ENERGY']:
        struc = dict['struc']
        weight = dict['weight']
        reaction = dict['reaction']
        divider = dict['divider']
        litVal = dict['litVal']
        reaxVal =0
        for j in range(1,len(reaction),2): #calculate reaxVal
            [subStruc,subDivider] = reaction[j].split('/')
            subDivider = int(subDivider)
            #print j,len(reaction)
            #print reaction[j],'test'
            if reaction[j-1] == '+':
                reaxVal += reaxDict[subStruc]['aEnergy']/subDivider
            elif reaction[j-1] == '-':
                reaxVal -= reaxDict[subStruc]['aEnergy']/subDivider
        #print reaction,reaxVal,litVal
        err = ((reaxVal-litVal)/weight)**2.0
        energyErr.append([weight,reaction,reaxVal,litVal,err])
    
    f.write('ENERGY\n')
    f.write('reaction\treaxVal\tlitVal\tweight\terror\ttotalError\n')
#    energyErr = sorted(energyErr,key = sortFun) #sorting
    for i in range(len(energyErr)):
        totalErr += energyErr[i][-1]
        energyErr[i].append(totalErr)
        line = energyErr[i]
            
        reaction = ' '.join(line[1])
        #reaction = reaction.replace('+ ','+ ').replace('- ','- ')
        f.write('%s\t%f\t%f\t%f\t%f\t%f\n' %(reaction,line[2],line[3],line[0],line[4],line[5]))

    #errDict = {'charge':chargeErr,'bond':bondErr,'angle':angleErr,'energy':energyErr}

    f.close()
    return totalErr
######

def readff(secID1,secID2,secID3):
    if secID1 == 1:
        startl=2
        startc=1
        nocol=1
        norow=1
    if secID1 == 2:
        startl=45 #start line
        startc=3 #start position
        nocol=8 #number of columns
        norow=4 #number of rows
    if secID1 ==3:
        startl=67
        startc=6
        nocol=8
        norow=2
    if secID1 == 4:
        startl=98
        startc=6
        nocol=6
        norow=1
    if secID1 == 5:
        startl=109
        startc=9
        nocol=7
        norow=1
    if secID1 == 6:
        startl=147
        startc=12
        nocol=7
        norow=1
    if secID1 == 7:
        startl=174
        startc=9
        nocol=4
        norow=1

    rowID = (secID2-1)*norow+startl+(secID3-1)//nocol #use python index, so no +1 here
    colID = (secID3-1)%nocol # same as above
    f = open('ffield','r')
    lines = f.readlines()
    f.close()
    line = lines[rowID]
    line = line[startc :]
    keywords = line.split()
    value = float(keywords[colID])
    return value
    
######
def changeff(secID1,secID2,secID3,value):
    if secID1 == 1:
        startl=2
        startc=1
        nocol=1
        norow=1
    if secID1 == 2:
        startl=45 #start line
        startc=3 #start position
        nocol=8 #number of columns
        norow=4 #number of rows
    if secID1 ==3:
        startl=67
        startc=6
        nocol=8
        norow=2
    if secID1 == 4:
        startl=98
        startc=6
        nocol=6
        norow=1
    if secID1 == 5:
        startl=109
        startc=9
        nocol=7
        norow=1
    if secID1 == 6:
        startl=147
        startc=12
        nocol=7
        norow=1
    if secID1 == 7:
        startl=174
        startc=9
        nocol=4
        norow=1

    rowID = (secID2-1)*norow+startl+(secID3-1)//nocol #use python index, so no +1 here
    colID = (secID3-1)%nocol # same as above
    f = open('ffield','r')
    lines = f.readlines()
    f.close()
    line = lines[rowID]
    startc += colID*9
    line = line[:startc]+'%9.4f' %(value)+line[startc+9:]
    lines[rowID] = line
#    print ''.join(lines)
    f = open('ffield','w')
    f.write(''.join(lines))
    f.close()

    


#####

def train():
    import shutil

    list = generate_inputs()
#############
   # print 'list'
   # for thing in list:
   #     print thing
   # print 'list end'
    reaxDict_old =reax_run(list)
    trainset_inDict = trainset_in()
    totError_old = calculate_error(reaxDict_old,trainset_inDict)
    shutil.copy('log.err','log.err.old')

    f = open('params','r')
    count = 0
    lines = f.readlines()
    for line in  lines:
      count += 1
      keywords = line.split()
      secID1 = int(keywords[0])
      secID2 = int(keywords[1])
      secID3 = int(keywords[2])
      delta = float(keywords[3])
      hiboundary = float(keywords[4])
      loboundary = float(keywords[5])

      print 'Parameter %d:(%d  %d  %d)' %(count,secID1,secID2,secID3)

      value = readff(secID1,secID2,secID3) #read original value
      if (value <= loboundary ) or (value  >= hiboundary):
          print '\tWarning: parameter out of boundary;%f not in [%f,%f] ' %(value,loboundary,hiboundary)
          print
          fdummy = open('warn_msg','a')
          line = '%d %d %d: Parameter out of boundary;%f not in [%f,%f]\n ' %(secID1,secID2,secID3,value,loboundary,hiboundary)
          fdummy.write(line)
          fdummy.close()

      else:
          trial1 = value-delta
          if trial1 <loboundary : trial1 = loboundary-0.05*abs(loboundary) #avoid same value
          #elif trial1 > hiboundary : trial1 = hiboundary+0.05*abs(hiboundary)
          changeff(secID1,secID2,secID3,trial1)
          reaxDict1 =reax_run(list)
          totError1 = calculate_error(reaxDict1,trainset_inDict)
          shutil.copy('log.err','log.err1')

          trial2 = value
          reaxDict2 = reaxDict_old
          totError2 = totError_old
          shutil.copy('log.err.old','log.err2')

          trial3 = value+delta
          if trial3 >hiboundary : trial3 = hiboundary+0.05*abs(hiboundary)
          #if trial3 <loboundary : trial3 = loboundary
          #elif trial3 > hiboundary : trial3 = hiboundary
          changeff(secID1,secID2,secID3,trial3)
          reaxDict3 =reax_run(list)
          totError3 = calculate_error(reaxDict3,trainset_inDict)
          shutil.copy('log.err','log.err3')

          a1 = totError1*(trial2-trial3)
          a2 = totError2*(trial3-trial1)
          a3 = totError3*(trial1-trial2)
          b1 = trial2-trial3
          b2 = trial3-trial1
          b3 = trial1-trial2
          c1 = trial2+trial3
          c2 = trial3+trial1
          c3 = trial1+trial2
          paraA = -(a1+a2+a3)/(b1*b2*b3)
          if paraA != 0.0: 
              value_opt = (a1*c1+a2*c2+a3*c3)/(a1+a2+a3)/2.0
          else: 
              value_opt = float('nan')
              print '\tWarning: possible dummy parameter'
              fdummy = open('warn_msg','a')
              line = '%d %d %d :Possible dummy parameter\n' %(secID1,secID2,secID3)
              fdummy.write(line)
              fdummy.close()

          reaxDict4 = {}
          totError4 = float('nan') #assign a big number
          trial4 = value_opt

          if paraA >0:
              if ((trial4 > loboundary )and (trial4 <hiboundary)):

                  changeff(secID1,secID2,secID3,trial4)
                  reaxDict4 =reax_run(list)
                  totError4 = calculate_error(reaxDict4,trainset_inDict)
                  shutil.copy('log.err','log.err4')

          if (totError4 <= totError1 and totError4<=totError2 and totError4<=totError3):
              reaxDict_old = reaxDict4
              totError_old = totError4
              shutil.copy('log.err4','log.err.old')
              newValue = trial4
          elif (totError2 <= totError1 and totError2 <= totError3):
              reaxDict_old = reaxDict2
              totError_old = totError2
              shutil.copy('log.err2','log.err.old')
              newValue = trial2
          elif (totError1 <= totError3):
              reaxDict_old = reaxDict1
              totError_old = totError1
              shutil.copy('log.err1','log.err.old')
              newValue = trial1
          else:
              reaxDict_old = reaxDict3
              totError_old = totError3
              shutil.copy('log.err3','log.err.old')
              newValue = trial3
          logFilename = 'log'+str(count)
          shutil.copy('log.err.old',logFilename)
          changeff(secID1,secID2,secID3,newValue)
          ffieldFile = 'ffield.run'+str(count)
          print '\trun1\trun2\trun3\trun4\tparaA'
          print '\t%f\t%f\t%f\t%f\t%f' %(trial1,trial2,trial3,trial4,paraA)
          print '\t%f\t%f\t%f\t%f' %(totError1,totError2,totError3,totError4)
          print '\tupdated value %f, backup as %s' %(newValue,ffieldFile)
          shutil.copy('ffield',ffieldFile)
    calculate_error(reaxDict_old,trainset_inDict)

def main():
    import shutil
    clean()

    #rewind warn_msg file
    f = open('warn_msg','w')
    f.close()

    args = sys.argv[1:]
    if not args:
        print 'using \'ffield\' file as starting point'
        train()
    elif len(args) == 1:
        print 'using \'' + args[0]+ '\'file as a starting point'
        inputFfield = str(args[0])
        shutil.copy(inputFfield,'ffield')
        train()
    
        
######
if __name__ == '__main__':
#    list=generate_inputs()
#    reaxDict = reax_run(list)
#    trainset_inDict = trainset_in()
#    totError = calculate_error(reaxDict,trainset_inDict)
#    readff(2,1,1)
#    readff(3,1,10)
#    readff(3,3,12)
#    changeff(2,3,8,3.5)
#    changeff(2,2,9,3.5)
#    print readff(2,4,7)
#    changeff(2,4,7,3.5)
    main()


    
