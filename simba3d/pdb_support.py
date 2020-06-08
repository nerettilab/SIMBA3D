
def make_pdb_file(outputfilepath,X):
    '''
    taks a 3xn array and converts it into pdb format
    '''
    smallest=abs(X).min()
    largest=abs(X).max()
    X*=100.0/largest
    with open(outputfilepath, 'w') as result:
        result.write("Structure Inference from Multiscale BAyes in 3 Dimensions (SIMBA3D)\n")
        for ii in range(len(X[0])):
            result.write('ATOM  ')
            result.write('{: 5d}'.format(ii+1))
            result.write('   CA MET A'+str(ii+1).ljust(8))
            result.write('{: 8.3f}'.format(X[0,ii]))
            result.write('{: 8.3f}'.format(X[1,ii]))
            result.write('{: 8.3f}'.format(X[2,ii]))
            result.write('  0.20 10.00\n')
        for ii in range(len(X[0])-1):
            result.write('CONECT')
            result.write('{: 5d}'.format(ii+1))
            result.write('{: 5d}'.format(ii+2)+'\n')