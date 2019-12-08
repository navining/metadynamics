def output(fileName, string):
    print(string.strip('\n'))
    f = open(fileName,'a')
    f.writelines(string)
    f.close()

def write_xyz(fileName, R):
    f = open('./xyz/'+fileName,'w')
    f.writelines(str(len(R))+'\n')
    f.writelines('Atom Positions\n')
    for atom in R:
        f.writelines('Ar\t'+str(atom[0])+'\t'+str(atom[1])+'\t'+str(atom[2])+'\n')
    f.close()

def write_R(fileName, R):
    f = open(fileName+'_R','w')
    for atom in R:
        f.writelines(str(atom[0])+', '+str(atom[1])+', '+str(atom[2])+'\n')
    f.close()

def write_V(fileName, V):
    f = open(fileName+'_V','w')
    for atom in V:
        f.writelines(str(atom[0])+', '+str(atom[1])+', '+str(atom[2])+'\n')
    f.close()
