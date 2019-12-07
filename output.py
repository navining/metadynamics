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



def write_xyz_(fileName, R):
    f = open('./xyz/'+fileName,'a+')
    f.writelines(str(len(R))+'\n')
    f.writelines('Atom Positions\n')
    for atom in R:
        f.writelines('Ar\t'+str(atom[0])+'\t'+str(atom[1])+'\t'+str(atom[2])+'\n')
    f.close()
    
