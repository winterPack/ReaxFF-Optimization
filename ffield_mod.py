__metaclass__ = type

class ffield:
    def __init__(self):
        self.paraNameSec1 = 'boc1,boc2,coa2,trip4,trip3,kc2,ovun6,trip2,ovun7,ovun8,trip1,swa,swb,n.u.,val6,lp1,val9,val10,n.u.,pen2,pen3,pen4,n.u.,tor2,tor3,tor4,n.u.,cot2,vdw1,cutoff,coa4,ovun4,ovun3,val8,n.u.,n.u.,n.u.,n.u.(gp37),coa3'.split(',')
        self.paraNameSec2 = 'ro_sigma,valency,mass,r_vdw,Dij,gammaEEM,ro_pi,val_e,alpha,gammaVdw,val_boc,ovun5,n.u.,chi,eta,p_hbond,ro_pipi,lp2,n.u.,boc4,boc3,boc5,n.u.,n.u.,ovun2,val3,n.u.,val_val,val5,rcore2,ecore2,acore2'.split(',')
        self.paraNameSec3 = 'De_sigma,De_pi,De_pipi,be1,bo5,13corr,bo6,ovun1,be2,bo3,bo4,n.u.,bo1,bo2,ovcorr,n.u.'.split(',')
        self.paraNameSec4 = 'Dij,Rdvw,alfa,ro_sigma,ro_pi,ro_pipi'.split(',')
        self.paraNameSec5 = 'theta0o,val1,val2,coa1,val7,pen1,val4'.split(',')
        self.paraNameSec6 = 'V1,V2,V3,tor1,cot1,n.u.,n.u.'.split(',')

    def say(self,secid):
        'output field value'
        paraNameSecs = [self.paraNameSec1,self.paraNameSec2,self.paraNameSec3,self.paraNameSec4,self.paraNameSec5,self.paraNameSec6]
        sec = 'sec'+str(secid[0])
        try:
            key = self.seq[secid[0]-1][secid[1]-1]            
            paraName = paraNameSecs[secid[0]-1][secid[2]-1]
            return ( (sec,key,paraName), self.ffDict[sec][key][paraName] )
        except IndexError:
            paraName = paraNameSecs[secid[0]-1][secid[2]-1]
            return ( (sec,paraName), self.ffDict[sec][paraName] )
            
    def mod(self,secid,newVal):
        'change ffield value'
        secinfo,par0 = self.say(secid)
        try:
            sec,key,paraName = secinfo
            self.ffDict[sec][key][paraName] = newVal
        except:
            sec,paraName = secinfo
            self.ffDict[sec][paraName] = newVal




    def set(self,filename):
        with open(filename,'r') as f:
            self.ffDict = {}
            self.seq = []

            #def readSec(f,secLabel,noRow):

            # read section 1
            # ln 1 to 41
            f.readline() #skip
            f.readline()
            buffer = []
            self.seq.append([])
            for ln in range(39):
                lnbuffer = f.readline()
                buffer.append(float(lnbuffer.split()[0]))
            #print buffer
            #dict(zip(paraNameSec1,buffer))
            #print zip(self.paraNameSec1,buffer)
            val = dict(zip(self.paraNameSec1,buffer))
            self.ffDict['sec1']=val;
            
            # read section 2
            tmpDict = {}
            self.seq.append([])
            lnbuffer = f.readline()
            NoLn = int(lnbuffer.split()[0])
            f.readline() #skip
            f.readline() #skip
            f.readline() #skip
            for i in range(NoLn):
                buffer = []
                for ln in range(4):
                    buffer = buffer+f.readline().split()
                key = buffer[0]
                self.seq[1].append(key)
                #print key
                val = dict(zip(self.paraNameSec2,map(float,buffer[1:])))
                tmpDict[key] = val
            self.ffDict['sec2'] = tmpDict

            # read section 3
            tmpDict = {}
            self.seq.append([])
            lnbuffer = f.readline()
            NoLn = int(lnbuffer.split()[0])
            f.readline() #skip
            for i in range(NoLn):
                buffer=[]
                for ln in range(2):
                    buffer = buffer + f.readline().split()
                key = tuple(buffer[0:2])
                self.seq[2].append(key)
                val = dict(zip(self.paraNameSec3,map(float,buffer[2:])))
                tmpDict[key] = val
            self.ffDict['sec3'] = tmpDict

            # read section 4
            tmpDict = {}
            self.seq.append([])
            lnbuffer = f.readline()
            NoLn = int(lnbuffer.split()[0])
            for i in range(NoLn) :
                buffer = f.readline().split()
                key = tuple(buffer[0:2])
                self.seq[3].append(key)
                val = dict(zip(self.paraNameSec4,map(float,buffer[2:])))
                tmpDict[key] = val
            self.ffDict['sec4'] = tmpDict

            # read section 5 
            tmpDict = {}
            self.seq.append([])
            lnbuffer = f.readline()
            NoLn = int(lnbuffer.split()[0])
            for i in range(NoLn):
                buffer = f.readline().split()
                key = tuple(buffer[0:3])
                self.seq[4].append(key)
                val = dict(zip(self.paraNameSec5,map(float,buffer[3:])))
                tmpDict[key] = val
            self.ffDict['sec5'] = tmpDict

            # read section 6
            tmpDict = {}
            self.seq.append([])
            lnbuffer = f.readline()
            NoLn = int(lnbuffer.split()[0])
            for i in range(NoLn):
                buffer = f.readline().split()
                key = tuple(buffer[0:4])
                self.seq[5].append(key)
                val = dict(zip(self.paraNameSec6,map(float,buffer[4:])))
                tmpDict[key] = val
            self.ffDict['sec6'] = tmpDict

            # read rest
            self.ffDict['rest'] = f.read()

            # set n.u. to 0
            for sec,secDict in self.ffDict.items():
                try:
                    for key in secDict.keys():
                        self.ffDict[sec][key]['n.u.'] = 0.0
                except AttributeError:
                    pass
                except TypeError:
                    self.ffDict[sec]['n.u.'] = 0.0







            

    def output(self):
        'sort and output force field to screen; obselete! use output2 instead'
        formatSec1 = '%10.4f !%s'
        formatSec2 = ' %-2s'+'%9.4f'*8 + '\n' + ' '*3 + '%9.4f'*8 + '\n   ' + '%9.4f'*8 + '\n   ' + '%9.4f'*8
        formatSec3 = '%3s'*2 + '%9.4f'*8 + '\n' + ' '*6+ '%9.4f'*8 
        formatSec4 = '%3s'*2 + '%9.4f'*6
        formatSec5 = '%3s'*3 + '%9.4f'*7
        formatSec6 = '%3s'*4 + '%9.4f'*7

        # print sec1
        print 'Reactive MD-force field\n 39        ! Number of general parameters'
        for i in range(39):
            key = self.paraNameSec1[i]
            print formatSec1 %(self.ffDict['sec1'][key],key)

        # print sec2
        print '%3d!' %len(self.ffDict['sec2']) + '%9s'*8 %tuple(self.paraNameSec2[0:8]) + '\n' \
        + ' '*4 + '%9s'*8 %tuple(self.paraNameSec2[8:16]) + '\n' \
        + ' '*4 + '%9s'*8 %tuple(self.paraNameSec2[16:24]) + '\n' \
        + ' '*4 + '%9s'*8 %tuple(self.paraNameSec2[24:32]) 
        for i in self.seq[1]:
            outlist = [i]
            for key in self.paraNameSec2:
                outlist.append(self.ffDict['sec2'][i][key])

            print formatSec2 % tuple(outlist)

        # print sec3
        print '%3d! ' %len(self.ffDict['sec3']) + '%9s'*8 %tuple(self.paraNameSec3[0:8]) + '\n' \
        + ' '*4 + '%9s'*8 %tuple(self.paraNameSec3[8:16])
        #print self.ffDict['sec3']
        for key in sorted(self.ffDict['sec3'].keys()):
            outlist = list(key)
            for key2 in self.paraNameSec3:
                val2 = self.ffDict['sec3'][key][key2]
                outlist.append(val2)
            print formatSec3 %tuple(outlist)

        # print sec4
        print '%3d   !' %len(self.ffDict['sec4'])+ '%9s'*6 %tuple(self.paraNameSec4)
        for key in sorted(self.ffDict['sec4'].keys()):
            outlist = list(key)
            for key2 in self.paraNameSec4:
                val2 = self.ffDict['sec4'][key][key2]
                outlist.append(val2)
            #print outlist
            print formatSec4 %tuple(outlist)

        # print sec5
       
        print '%3d   !' %len(self.ffDict['sec5'])+ '%9s'*7 %tuple(self.paraNameSec5)
        for key in sorted(self.ffDict['sec5'].keys()):
            outlist = list(key)
            for key2 in self.paraNameSec5:
                val2 = self.ffDict['sec5'][key][key2]
                outlist.append(val2)
            print formatSec5 %tuple(outlist)

        # print sec6
        print '%3d     !' %len(self.ffDict['sec6']) + '%9s'*7 %tuple(self.paraNameSec6)
        for key in sorted(self.ffDict['sec6'].keys()):
            outlist = list(key)
            for key2 in self.paraNameSec6:
                val2 = self.ffDict['sec6'][key][key2]
                outlist.append(val2)
            print formatSec6 %tuple(outlist)

        # print rest
        print self.ffDict['rest']

    def output2file(self,filename):
        'sort and output force field to file'
        content = ''
        formatSec1 = '%10.4f !%s\n'
        formatSec2 = ' %-2s'+'%9.4f'*8 + '\n' + ' '*3 + '%9.4f'*8 + '\n   ' + '%9.4f'*8 + '\n   ' + '%9.4f'*8 + '\n'
        formatSec3 = '%3s'*2 + '%9.4f'*8 + '\n' + ' '*6+ '%9.4f'*8 + '\n'
        formatSec4 = '%3s'*2 + '%9.4f'*6 + '\n'
        formatSec5 = '%3s'*3 + '%9.4f'*7 + '\n'
        formatSec6 = '%3s'*4 + '%9.4f'*7 + '\n'

        # print sec1
        content = content + 'Reactive MD-force field\n 39        ! Number of general parameters\n'
        for i in range(39):
            key = self.paraNameSec1[i]
            content = content + formatSec1 %(self.ffDict['sec1'][key],key)

        # print sec2
        content = content + '%2d !' %len(self.ffDict['sec2']) + '%9s'*8 %tuple(self.paraNameSec2[0:8]) + '\n' \
        + ' '*4 + '%9s'*8 %tuple(self.paraNameSec2[8:16]) + '\n' \
        + ' '*4 + '%9s'*8 %tuple(self.paraNameSec2[16:24]) + '\n' \
        + ' '*4 + '%9s'*8 %tuple(self.paraNameSec2[24:32]) + '\n'
        for i in self.seq[1]:
            outlist = [i]
            for key in self.paraNameSec2:
                outlist.append(self.ffDict['sec2'][i][key])

            content = content + formatSec2 % tuple(outlist)

        # print sec3
        content = content + '%2d ! ' %len(self.ffDict['sec3']) + '%9s'*8 %tuple(self.paraNameSec3[0:8]) + '\n' \
        + ' '*4 + '%9s'*8 %tuple(self.paraNameSec3[8:16]) + '\n'
        #print self.ffDict['sec3']
        for key in self.seq[2]:
            outlist = list(key)
            for key2 in self.paraNameSec3:
                val2 = self.ffDict['sec3'][key][key2]
                outlist.append(val2)
            content = content + formatSec3 %tuple(outlist)

        # print sec4
        content = content + '%3d   !' %len(self.ffDict['sec4'])+ '%9s'*6 %tuple(self.paraNameSec4) + '\n'
        for key in self.seq[3]:
            outlist = list(key)
            for key2 in self.paraNameSec4:
                val2 = self.ffDict['sec4'][key][key2]
                outlist.append(val2)
            #print outlist
            content = content + formatSec4 %tuple(outlist)

        # print sec5
       
        content = content + '%3d   !' %len(self.ffDict['sec5'])+ '%9s'*7 %tuple(self.paraNameSec5) + '\n'
        for key in self.seq[4]:
            outlist = list(key)
            for key2 in self.paraNameSec5:
                val2 = self.ffDict['sec5'][key][key2]
                outlist.append(val2)
            content = content + formatSec5 %tuple(outlist)

        # print sec6
        content = content + '%3d     !' %len(self.ffDict['sec6']) + '%9s'*7 %tuple(self.paraNameSec6) + '\n'
        for key in self.seq[5]:
            outlist = list(key)
            for key2 in self.paraNameSec6:
                val2 = self.ffDict['sec6'][key][key2]
                outlist.append(val2)
            content = content + formatSec6 %tuple(outlist)

        # print rest
        content = content + self.ffDict['rest']
        with open(filename,'w') as f:
            f.write(content)


    def cmp(self,ff):
        'compare with a reference force field file, using supplied file as an ref'
        if not isinstance(ff,ffield):
            raise 'instance of ffield class needed'

        cmpDict = {}
        for secid in 'sec1,sec2,sec3,sec4,sec5,sec6'.split(','):
            cmpDict[secid]={}
            for key in ff.ffDict[secid].keys():
                try:
                    val = {}
                    for key2 in ff.ffDict[secid][key].keys():
                        val2Ref = ff.ffDict[secid][key][key2]
                        val2 = self.ffDict[secid][key][key2]
                        try:
                            valDif = (val2-val2Ref, (val2-val2Ref)/val2Ref*100 )
                        except ZeroDivisionError:
                            valDif = (val2-val2Ref, float('nan') )
                        val[key2] = valDif

                except AttributeError:
                    try:
                        val = ( self.ffDict[secid][key]-ff.ffDict[secid][key] , self.ffDict[secid][key]/ff.ffDict[secid][key]*100-100)
                    except ZeroDivisionError:
                        val = ( self.ffDict[secid][key]-ff.ffDict[secid][key], float('nan') )
                    #print self.ffDict[secid][key],ff.ffDict[secid][key],val
                except KeyError:
                    val = None
                cmpDict[secid][key] = val


        # print cmp result
        # print sec1
        print '--- Section 1 ---'
        print '%-19s' %'para Name',
        print '%12s%12s   %12s'  %('ref Value','difference','   percentile')

        count = 0
        for key in ff.paraNameSec1:
            count = count + 1
            print 'param %2d: %-9s' %(count,key) + '%12.4f' %ff.ffDict['sec1'][key] + '%+12.4f (%12.2f%%)' %cmpDict['sec1'][key]        

        # print sec2
        print '--- Section 2 ---'
        print '%-19s' %'para Name',
        print '|            %12s                  '*len(ff.ffDict['sec2']) %tuple(ff.seq)
        count = 0
        for key in ff.paraNameSec2:
            print 'param %2d: %-9s' %(count,key) ,
            for key2 in ff.seq[1]:
                print '|%12.4f' %ff.ffDict['sec2'][key2][key],
                try:
                    print '%+12.4f (%12.2f%%)' %cmpDict['sec2'][key2][key] ,
                except TypeError:
                    print '%12.4f (%12.2f%%)' %(float('nan') , float('nan')) ,
            print '|'

        # print sec3
        print '--- Section 3 ---'
        print '%-19s' %'para Name',
        for tup in sorted(ff.ffDict['sec3'].keys()):
            print '|             %3s    ,%3s                 ' %tup,
        print '|'

        count = 0
        for key in ff.paraNameSec3:
            count = count + 1
            print 'param %2d: %-9s' %(count,key) ,
            for key2 in self.seq[2]:
                print '|%12.4f' %ff.ffDict['sec3'][key2][key],
                try:
                    print '%+12.4f (%12.2f%%)' %cmpDict['sec3'][key2][key] ,
                except TypeError:
                    print '%12.4f (%12.2f%%)' %(float('nan') , float('nan')) ,
            print '|'

        # print sec4
        print '--- Section 4 ---'
        print '%-19s' %'para Name',
        for tup in sorted(ff.ffDict['sec4'].keys()):
            print '|             %3s    ,%3s                 ' %tup,
        print '|'

        count = 0
        for key in ff.paraNameSec4:
            count = count + 1
            print 'param %2d: %-9s' %(count,key) ,
            for key2 in self.seq[3]:
                print '|%12.4f' %ff.ffDict['sec4'][key2][key],
                try:
                    print '%+12.4f (%12.2f%%)' %cmpDict['sec4'][key2][key] ,
                except TypeError:
                    print '%12.4f (%12.2f%%)' %(float('nan') , float('nan')) ,
            print '|'

        # print sec5
        print '--- Section 5 ---'
        print '%-19s' %'para Name',
        for tup in sorted(ff.ffDict['sec5'].keys()):
            print '|          %3s  ,%3s  ,%3s                ' %tup,
        print '|'

        count = 0
        for key in ff.paraNameSec5:
            count = count + 1
            print 'param %2d: %-9s' %(count,key) ,
            for key2 in self.seq[4]:
                print '|%12.4f' %ff.ffDict['sec5'][key2][key],
                try:
                    print '%+12.4f (%12.2f%%)' %cmpDict['sec5'][key2][key] ,
                except TypeError:
                    print '%12.4f (%12.2f%%)' %(float('nan') , float('nan')) ,
            print '|'

        # print sec6
        print '--- Section 6 ---'
        print '%-19s' %'para Name',
        for tup in sorted(ff.ffDict['sec6'].keys()):
            print '|       %3s  ,%3s  ,%3s  ,%3s             ' %tup,
        print '|'

        count = 0
        for key in ff.paraNameSec6:
            count = count + 1
            print 'param %2d: %-9s' %(count,key) ,
            for key2 in self.seq[5]:
                print '|%12.4f' %ff.ffDict['sec6'][key2][key],
                try:
                    print '%+12.4f (%12.2f%%)' %cmpDict['sec6'][key2][key] ,
                except TypeError:
                    print '%12.4f (%12.2f%%)' %(float('nan') , float('nan')) ,
            print '|'


class LAMMPS:
    'all works needed to be done by LAMMPS are here'
    def __init__(self):
        self.generate_inputs()

    def generate_inputs(self):
        'generate inputs for lammps run'
        fileIn1 = open('geo','r')
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
                    outDict = { 'massList':massList, 'elementList':elementList,'cName':cName}
                    content = '''
clear
log        log.%(cName)s
#echo            both

units           real

atom_style      charge
read_data       %(cName)s.str
%(massList)s

pair_style      reax/c NULL

pair_coeff      * * ffield.tmp %(elementList)s

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
#dump            1 all atom 1 dump.reaxc
# dump file 2
# dump file 3
#dump        5 all xyz 1 dump.xyz 
# dump file 4
dump        10 all custom 1 dump.custom.%(cName)s id element q x y z
dump_modify     10 element %(elementList)s

minimize        1e-6 1e-8 10000 100000''' % (outDict)
                    inFile.write(content)
                    inFile.close()

        self.structures = structures

    def runAll(self):
        #import os
        import subprocess
        from math import ceil
        noStructures = len(self.structures)
        NPROC = 8
        stepSize = int( ceil ( noStructures/float(NPROC) ))
        cmdLine = ''
        for i in range(0,noStructures,stepSize):
            try:
                #print i+stepSize
                #print self.structures[i:i+stepSize]
                cmdLine = cmdLine + ' cat ' + ' '.join( map(lambda x: 'in.'+ x , self.structures[i:i+stepSize])) + ' |./lmp_openmpi >/dev/null &'
            except:
                raise

        cmdLine = cmdLine + ' wait'
        #print ' '.join(self.structures),len(self.structures), stepSize
        #print cmdLine

        try:
            subprocess.check_call(cmdLine,shell=True)
        except subprocess.CalledProcessError:
            print '!!! LAMMPS run error !!!' 
            raise

    def clean(self):
        import os
        cmdLine = [' dump.custom.%s log.%s' %(x,x) for x in self.structures]
        cmdLine = 'rm -rf ' + ' '.join(cmdLine)
        os.system(cmdLine)

    def readDumpLog(self):
        import re

        self.minDict = {} #info about minimized structures
        for cName in self.structures:
            dumpFile = 'dump.custom.%s' %cName
            #print dumpFile,'dumpFile name'
            tmpDict = {}
            f = open(dumpFile,'r') #dump file
            content = f.readlines()
            f.close()
            noa = int(content[3])
            nol = len(content)
            for line in content[-noa:]:
                keys = line.split()
                tmpDict2 = {}
                tmpDict2['coordination'] = map(float,keys[-3:])
                tmpDict2['q'] = float(keys[2])
                tmpDict2['type'] = keys[1]
                tmpDict[keys[0]] = tmpDict2
            logFile = 'log.' + cName
            f = open(logFile,'r') #log file
            content = f.read()
            f.close()
            eng = float(re.findall(r'\n\s*\d+\s+0\s+([-0-9.]+)\s+.*\nLoop\stime\sof\s',content)[0])
            #print eng, cName
            tmpDict['energy'] = eng
            self.minDict[cName] = tmpDict

    def sayCharge(self,cName,index):
        'output is charge info'
        return self.minDict[cName][index]['q']

    def sayGeometry(self,cName,info):
        'output distance or angle'
        import math
        dotProd = lambda vec1,vec2 : sum([vec1[x]*vec2[x] for x in range(3)])
        vec = lambda cName,index1,index2: [self.minDict[cName][index2]['coordination'][x]- self.minDict[cName][index1]['coordination'][x] for x in range(3)]
        if len(info) == 2:
            index1 = info[0]
            index2 = info[1]
            return dotProd(vec(cName,index1,index2),vec(cName,index1,index2))**0.5
        elif len(info) == 3:
            index1 = info[0]
            index2 = info[1]
            index3 = info[2]
            vec1 = vec(cName,index2,index1)
            vec2 = vec(cName,index2,index3)
            cosVal = dotProd(vec1,vec2) / (dotProd(vec1,vec1)*dotProd(vec2,vec2))**0.5 
            cosVal = min(cosVal,1)
            cosVal = max(cosVal,-1)
            return math.degrees(math.acos(cosVal))

    def sayEnergy(self,cName):
        'output is energy'
        return self.minDict[cName]['energy']




class trainset:
    def __init__(self):
        self.trainDict = {'CHARGE':[], 'GEOMETRY':[], 'ENERGY':[]}
        self.outContent = {} # store output info
        #print self.trainDict
        self.readfile()
    def readfile(self):
        'read trainset.in file info'
        with open('trainset.in','r') as f:
            for line in f:
                if (line.strip()== '') or (line.strip()[0] == '#'): continue
                words = line.split()
                if words[0][0:3] == 'END': continue

                if len(words) == 1:
                    key = words[0]
                else:
                    tmpDict = {'litVal':float(words[-1]),'reaxVal':None, 'error':None}
                    try:
                        tmpDict.update({'cName':words[0], 'weight':float(words[1]),'info':words[2:-1]})
                    except ValueError:
                        tmpDict.update({'weight':float(words[0]),'info':words[1:-1]})
                    self.trainDict[key].append(tmpDict)
        # set torlerance
        # charge
        for n in range(len(self.trainDict['CHARGE'])):
            self.trainDict['CHARGE'][n].update({'torlerance':10.0})
        # geo
        for n in range(len(self.trainDict['GEOMETRY'])):
            info = self.trainDict['GEOMETRY'][n]['info']
            if len(info) == 2:
                torlerance = 0.1 # unit: A
            if len(info) == 3:
                torlerance = 5 # unit: degree
            self.trainDict['GEOMETRY'][n].update({'torlerance':torlerance})
        # energy
        for n in range(len(self.trainDict['ENERGY'])):
            self.trainDict['ENERGY'][n].update({'torlerance':30.0})




    def checkError(self,myLAMMPS):
        'check error against literature value'
        # update charge vals
        for n in range(len(self.trainDict['CHARGE'])):
            line = self.trainDict['CHARGE'][n]
            cName = line['cName']
            index = line['info'][0]
            reaxVal = myLAMMPS.sayCharge(cName,index)
            self.trainDict['CHARGE'][n]['reaxVal'] = reaxVal
            #print reaxVal
        # update geo vals
        for n in range(len(self.trainDict['GEOMETRY'])):
            line = self.trainDict['GEOMETRY'][n]
            cName = line['cName']
            info = line['info']
            reaxVal = myLAMMPS.sayGeometry(cName,info)
            self.trainDict['GEOMETRY'][n]['reaxVal'] = reaxVal
            #print reaxVal
        # update energy vals
        for n in range(len(self.trainDict['ENERGY'])):
            info = self.trainDict['ENERGY'][n]['info']
            reaxVal = 0
            for k in range(0,len(info),2):
                sign = float( info[k] + '1')
                words = info[k+1].split('/')
                cName = words[0]
                quantity = float(words[1])
                reaxVal =  reaxVal + sign * myLAMMPS.sayEnergy(cName)/quantity
                self.trainDict['ENERGY'][n]['reaxVal'] = reaxVal

        # calculate error
        for key,value in self.trainDict.items():
            for n in range(len(value)):
                self.trainDict[key][n]['error'] = self.trainDict[key][n]['reaxVal'] - self.trainDict[key][n]['litVal']

    def fatError(self,myLAMMPS):
        'use a filter function to sum error'
        self.checkError(myLAMMPS)
        import copy
        heaviside = lambda x: 1.0 if x >0.0 else -0.0
        fatCalc = lambda x,torlerance: heaviside(abs(x) - torlerance) * (abs(x) - torlerance)
        dictCalc = lambda dict: (fatCalc(dict['error'],dict['torlerance']) / dict['weight'])**2.0
        for key in self.trainDict.keys():
            for n in range(len(self.trainDict[key])):
                line = self.trainDict[key][n]
                fatErrorVal = dictCalc(line)
                self.trainDict[key][n].update({'fatError':fatErrorVal})
        # generate output info
        content = ''
        colHeader = '%-40s' + '%-16s'*7 %('weight','litVal','reaxVal','error','torlerance','fatError(weighted)','sumError(weighted)') + '\n'
        errSum = 0
        for key in ('CHARGE', 'GEOMETRY', 'ENERGY'):
            content = content + colHeader %key
            for line in self.trainDict[key]:
                errSum = errSum + line['fatError']
                outDict = copy.deepcopy(line)
                header = ''
                try:
                    cName = outDict['cName']
                    header = header + '%-16s' %cName
                    for n in outDict['info']:
                        header = header + '%2s/%-2s ' %(n,myLAMMPS.minDict[cName][n]['type'])
                except:
                    header = ' '.join(outDict['info'])
                outDict['header'] = header
                outDict['errSum'] = errSum

                content = content + '%(header)-40s%(weight)-16.4f%(litVal)-16.4f%(reaxVal)-16.4f%(error)-16.4f%(torlerance)-16.4f%(fatError)-16.4f%(errSum)-16.4f\n' %outDict
        self.outContent.update({'fatError':content})
        return errSum

        


                    

class params:
    def __init__(self):
        self.params = []
        with open('params','r') as f:
            for line in f:
                if ( line.strip() != '') and (line.strip()[0] != '#'):
                    words = line.split()
                    self.params.append( [map(int,words[0:3]), map(float,words[3:6])] )
        self.nl = len(self.params)



                


class reaxFit:
    def __init__(self):
        self.myLAMMPS = LAMMPS()
        self.myFfield = ffield()
        self.myTrainset = trainset()





if __name__ == "__main__":
    import sys
    buffer = sys.argv[1:]
    if len(buffer) == 1:
        ffname = buffer[0]
        ff1 = ffield()
        ff1.set(ffname)
        ff1.output()
    elif len(buffer) == 2:
        ff1 = ffield()
        ff2 = ffield()
        ffname1,ffname2 = buffer
        ff1.set(ffname1)
        ff2.set(ffname2)
        print 'Check %s use %s as a refference' %(ffname2,ffname1)
        ff2.cmp(ff1)
    elif len(buffer) == 0:
        import copy

        def drange(start,finish,delta):
            'range() function in floats; return a list'
            assert (finish-start)*delta > 0
            res = []
            i = start
            sign = round(delta/abs(delta))
            while i * sign < finish * sign:
                res.append(i)
                i = i + delta
            return res

        def splot(xs,ys):
            'simple plot'
            ny=5.0
            xmin = xs[0]
            xmax = xs[-1]
            ymin = min(ys)
            ymax = max(ys)
            delta = (ymax-ymin)/ny
            ylabel = drange(ymin,ymax+0.000001,delta)
            print '\n'
            print '%-16s' %'***EzPlot' +'%-16.4f'*len(xs) %tuple(xs) 
            for nyl in range(len(ylabel)):
                line = '%-16.4f' %ylabel[nyl]
                for y in ys:
                    if round( (y-ymin)/delta ) == nyl:
                        line = line + '*' + ' '*15 
                    else:
                        line = line + ' '*16
                print line
        def matlibPlt(xs,ys,text,filename):
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot(xs,ys,'bo-')
            plt.title(text)
            plt.xlabel('param')
            plt.ylabel('error')
            plt.grid(True)
            plt.savefig(filename)


        # initialization
        myParams = params()
        myReaxFit = reaxFit()
        myReaxFit.myFfield.set('ffield')
        myReaxFit.myFfield.output2file('ffield.tmp')
        myReaxFit.myFfield.output2file('ffield.0')
        myReaxFit.myLAMMPS.generate_inputs()
        myReaxFit.myLAMMPS.runAll()
        myReaxFit.myLAMMPS.readDumpLog()
        myReaxFit.myTrainset.fatError(myReaxFit.myLAMMPS)
        with open('log.err.0','w') as f:
            f.write(myReaxFit.myTrainset.outContent['fatError'])

        count = 0
        for secid,parConfig in myParams.params:
            count = count + 1
            testReaxFit = copy.deepcopy(myReaxFit)
            delta,high,low = parConfig
            secinfo,par0 = testReaxFit.myFfield.say(secid)
            print '\n'*2
            print 'Section ID:' + '%4d'*3 %tuple(secid) + '    Value:%16.4f' %par0 + '    delta:%16.4f    Range:    [%16.4f,%16.4f]' %(delta,low,high)
            print '***',secinfo
            parTrial = drange(par0,low,-delta)[1:]
            parTrial.reverse()
            parTrial = parTrial + drange(par0,high,delta)
            errList = []
            print '%-16s' %'parTrial:' + '%-16.4f'*len(parTrial) %tuple(parTrial)
            sys.stdout.write('%-16s' %'error:')
            for par in parTrial:
                testReaxFit.myFfield.mod(secid,par)
                testReaxFit.myFfield.output2file('ffield.tmp')
                testReaxFit.myLAMMPS.runAll()
                testReaxFit.myLAMMPS.readDumpLog()
                #print testReaxFit.myTrainset.trainDict
                err = testReaxFit.myTrainset.fatError(testReaxFit.myLAMMPS)
                errList.append(err)
                # print '%-16s' %'error:' +  '%-16.4f'*len(errList) %tuple(errList) + '\r',
                sys.stdout.write('%-16.4f' %err)
                sys.stdout.flush()
            print
            text = '%s %s' %(str(secid),str(secinfo))
            figFile = 'plot%d.png' %count
            #splot(parTrial,errList)
            matlibPlt(parTrial,errList,text,figFile) # plot
            minIndex = [i for i,x in enumerate(errList) if x == min(errList)]
            minErr = [errList[i] for i in minIndex ]
            newIndex = minIndex[int(len(minIndex)/2.0+0.0001) - 1]
            newPar = parTrial[newIndex]
            newError = errList[newIndex]
            if newPar != par0:
                myReaxFit.myFfield.mod(secid,newPar)
                print '\n***' + 'Value updated!    New Value: %-16.4f    New Error: %-16.4f' %(newPar,newError)
                myReaxFit.myFfield.output2file('ffield.tmp')
            else:
                print '\n***' + 'No updates!'
            myReaxFit.myFfield.output2file('ffield.'+str(count))
            myReaxFit.myLAMMPS.runAll()
            myReaxFit.myLAMMPS.readDumpLog()
            myReaxFit.myTrainset.fatError(myReaxFit.myLAMMPS)
            with open('log.err.'+str(count),'w') as f:
                f.write('log file for ffield.%d\n' %count)
                f.write(myReaxFit.myTrainset.outContent['fatError'])








                
            # for par in parTrial:

        # myReaxFit.myLAMMPS.runAll()
        # myReaxFit.myLAMMPS.readDumpLog()
        # myReaxFit.myTrainset.checkError(myReaxFit.myLAMMPS)
        # myReaxFit.myTrainset.fatError(myReaxFit.myLAMMPS)


        # for line in myParams.params:
        #     secid = line[0]
        #     info = line[1]
        #     print secid

            #print myReaxFit.myFfield.say(secid)
            

        # myLAMMPS = LAMMPS()
        # myLAMMPS.runAll()
        # myLAMMPS.readDumpLog()
        # #myLAMMPS.clean()
        # myTrainset = trainset()
        # myTrainset.checkError(myLAMMPS)
        # myFfield = ffield()
        # myFfield.set('ffield')
        # print myFfield.say([2,1,6])
        # myff1 = ffield()
        # myff1.set('ffield')
        # print myff1.say([2,1,1])
        # myff2 = copy.deepcopy(myff1)
        # myff2.mod([2,1,1],99)
        # print myff1.say([2,1,1])
        # print myff2.say([2,1,1])
        # myff2.output2file('test')
        # #print myff2.ffDict['sec1']
        # # myff2.output()
        # print myff2.say([2,1,1])
        # #print myff2.ffDict






