from train import *

def read_params():
    filename = 'params'
    params = []
    with open(filename,'r') as f:
        for line in f:
            keywords = line.split()
            secID1 = int(keywords[0])
            secID2 = int(keywords[1])
            secID3 = int(keywords[2])
#            info = line[line.index('!'):]
            info = ' ! auto_gen\n'
            params.append([secID1,secID2,secID3,info])
    return params


def generate_params(params,scale):
    filename = 'params.new'
    with open(filename,'w') as f:
        for list in params:
            secID1 = list[0]
            secID2 = list[1]
            secID3 = list[2]
            info = list[3]
            value = readff(secID1,secID2,secID3)
            delta = abs(value*scale)
            if delta <= 0.0001: delta = 0.0001
            dbound = 2.0*delta
            loboundary = value-dbound
            hiboundary = value+dbound
            #line = '%d\t%d\t%d\t%f\t%f\t%f\t%s' %(secID1,secID2,secID3,delta,hiboundary,loboundary,info)
            line = '%3d%3d%3d%16.4f%16.4f%16.4f %s' %(secID1,secID2,secID3,delta,hiboundary,loboundary,info)
            f.write(line)

def auto_train():
    import sys
    import shutil
    args = sys.argv[1:]
    scale = float(args[0])
    #rewind warn_msg file
    print 'auto train start from scale = %f' %scale
    for i in range(3):
        print 'run %d:' %(i)
        params = read_params()
        generate_params(params,scale)
        shutil.copy('params.new','params')
        train()    




def main():
    f = open('warn_msg','w')
    f.close()
    import sys
    args = sys.argv[1:]
    scale = float(args[0])
    print 'using scale = %f to generate params.new file' %(scale)
    params = read_params()
    generate_params(params,scale)
   # auto_train()

if __name__ == '__main__':
    main()
