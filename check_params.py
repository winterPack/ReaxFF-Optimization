import train

f = open('params','r')
lines = f.readlines()

for line in lines:
    keywords = line.split()
    secID1 = int(keywords[0])
    secID2 = int(keywords[1])
    secID3 = int(keywords[2])
    value = train.readff(secID1,secID2,secID3)
    hibound = float(keywords[4])
    lobound = float(keywords[5])
    if (value > lobound) and (value <hibound):
        check = 'ok'
    else: 
        check = 'no'


    print check, secID1,secID2,secID3,value,hibound,lobound,line[39:],

    
