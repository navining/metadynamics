def output(fileName, string):
    print(string.strip('\n'))
    f = open(fileName,'a')
    f.writelines(string)
    f.close()