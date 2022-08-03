from __future__ import print_function
import numpy
import scipy
import pdb

''' this function takes one electron and two electron integral files and makes new integral files in which the basis functions are eigenfunctions of 
D_{inf}h point group. '''
'''
Lz0 = [0, 1, 4, 5, 6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 2, 3, 38, 39, 40, 41, 42, 44, 45, 47, 48, 49, 50, 51]
Lz1a = [[20, 21, 22, 23, 24, 25, 26, 27], [28, 29, 30, 31, 32, 33, 34, 35]]
Lz1b = [[52, 53, 54, 55, 56, 57, 58, 59], [60, 61, 62, 63, 64, 65, 66, 67]]
Lz2 = [[9, 13, 36, 37], [18, 19, 43, 46]]
Larray = [Lz0, Lz1a, Lz1b, Lz2];
mLz = [0,1,1,2]
'''
''' pvdz'''
norbs = 68
Lz0 = [0, 2, 6, 7, 8, 10, 11, 13, 14, 16, 34, 36, 37, 38, 40, 41, 42, 44, 45, 47]
Lz1a = [[19, 17, 18, 3, 20, 21, 22], [26, 24, 25, 4, 27, 28, 29]]
Lz1b = [[48, 49, 50, 51, 52, 53, 55], [56, 57, 58, 59, 60, 61, 63]]
Lz2a = [[1, 12, 9, 15], [5, 32, 31, 33]]
Lz2b = [[35, 43, 46, 39], [64, 66, 67, 65]]
Lz3a = [[23], [30]]
Lz3b = [[54], [62]]
#Lz4a = [[15, 26], [43, 48]]
#Lz4b = [[121, 131], [101, 105]]
Larray = [Lz0, Lz1a, Lz1b, Lz2a, Lz2b, Lz3a, Lz3b]#, Lz4a, Lz4b];
mLz = [0,1,1,2,2,3,3]#,4, 4]

#negatives = [84,85,82,  125,126,128,129,131,132,133,135,136,138,139,140,142,143,144,  39,40,41,45,  86,90,91,  79]
negatives = []#[92,94,  161,165,167,168,   40,44,45,46,47,49,  108,  80,89,90,  177,  101]

screen_error = 1e-9

int2 = numpy.zeros(shape=(norbs, norbs, norbs, norbs));
int1 = numpy.zeros(shape=(norbs, norbs));

coeffs = numpy.zeros(shape=(norbs, norbs)).astype(complex);

index = 0;
m = 0;
lindex = 0
for Li in Larray:
   m = mLz[lindex]
   if type(Li[0]) == type(1):
       print("transforming orbs with L = 0")
       for i in Li:
           coeffs[index, i] = 1.0
           index =index+1
   else :
       for j in range(len(Li[0])):
           scale = 1.0
           if (Li[1][j] in negatives):
              scale = -1.0
              print("index with negatives ", index, Li[0][j], "+j", Li[1][j])
              print("index with negatives ", index+1, Li[0][j], "-j", Li[1][j])
           coeffs[index, Li[0][j] ] = 1.0/(2.0**0.5)
           #coeffs[index, Li[1][j] ] = -1.0/(2.0**0.5)
           coeffs[index, Li[1][j] ] = scale*-1.0j/(2.0**0.5)
           index = index+1

           coeffs[index, Li[0][j] ] = 1.0/(2.0**0.5)
           #coeffs[index, Li[1][j] ] = ((-1)**m)*1.0/(2.0**0.5)
           coeffs[index, Li[1][j] ] = scale*1.0j/(2.0**0.5)
           index = index+1
   lindex=lindex+1



print("reading one ints")
f = open("qcdmrg.int1","r")
lines = f.readlines()
n = int(lines[0].split()[0])
if (n != norbs) :
    exit(-1)
for line in lines[1:]:
    tokens = line.split()
    int1[int(tokens[0]),int(tokens[1])] = float(tokens[2])
    int1[int(tokens[1]),int(tokens[0])] = float(tokens[2])

print("one int transform")
#newint1 = numpy.dot(numpy.dot(coeffs,int1),coeffs.conj().transpose())
newintt = numpy.tensordot(coeffs.conj(), int1, axes=([1],[0]))
newint1 = numpy.tensordot(newintt, coeffs, axes=([1],[1]))
f.close()

f = open("new.int1","w")
f.write(str(norbs)+'\n')
for i in range(norbs):
  for j in range(norbs):
     if (abs(newint1[i,j]) >=screen_error):
        if (newint1[i,j].imag >= screen_error):
           print('problem in ', i, j, newint1[i,j])
        s = str(i)+'  '+str( j)+"  "+str( newint1[i,j].real)+'\n'
        f.write(s)
f.close()



print("read two ints")
#print numpy.dot(coeffs, coeffs.conj().transpose())
#now read integrals
f = open("qcdmrg.int2","r")
lines = f.readlines()
n = int(lines[0].split()[0])
if (n != norbs) :
    print("number of orbitals don't match", n, norbs)
    exit(-1)
for line in lines[1:]:
    tokens = line.split()
    int2[int(tokens[0]),int(tokens[1]),int(tokens[2]),int(tokens[3])] = float(tokens[4])
    int2[int(tokens[2]),int(tokens[1]),int(tokens[0]),int(tokens[3])] = float(tokens[4])
    int2[int(tokens[0]),int(tokens[3]),int(tokens[2]),int(tokens[1])] = float(tokens[4])
    int2[int(tokens[2]),int(tokens[3]),int(tokens[0]),int(tokens[1])] = float(tokens[4])

    int2[int(tokens[1]),int(tokens[0]),int(tokens[3]),int(tokens[2])] = float(tokens[4])
    int2[int(tokens[3]),int(tokens[0]),int(tokens[1]),int(tokens[2])] = float(tokens[4])
    int2[int(tokens[1]),int(tokens[2]),int(tokens[3]),int(tokens[0])] = float(tokens[4])
    int2[int(tokens[3]),int(tokens[2]),int(tokens[1]),int(tokens[0])] = float(tokens[4])

f.close()


newintt = numpy.zeros(shape=(norbs, norbs, norbs, norbs));
print("step1")

#numpy.einsum('ij,jklm->iklm', coeffs.conj(), int2, newintt)
newintt = numpy.tensordot(coeffs.conj(), int2, axes=([1],[1]))
print("step2")
#numpy.einsum('ij,kjlm->kilm', coeffs.conj(), newintt, int2)
newint2 = numpy.tensordot(coeffs.conj(), newintt, axes=([1],[1]))
print("step3")
#numpy.einsum('ij,kljm->klim', coeffs, int2, newintt)
newintt = numpy.tensordot(newint2, coeffs, axes=([2],[1]))
print("step4")
#numpy.einsum('ij,klmj->klmi', coeffs, newintt, int2)
newint2 = numpy.tensordot(newintt, coeffs, axes=([2],[1]))

print("writing ints")

f = open("new.int2","w")
f.write(str(norbs)+'\n')
for i in range(norbs):
  for j in range(norbs):
      for k in range(norbs):
          for l in range(norbs):
             if (abs(newint2[i,j,k,l]) >=screen_error ):
                if (newint2[i,j,k,l].imag >= screen_error):
                   print('problem in ', i, j, k,l,newint2[i,j,k,l])

                s = str(i)+'  '+str( j)+"  "+str( k)+"  "+str( l)+"  "+str( newint2[i,j, k, l].real)+'\n'
                f.write(s)
f.close()

'''
err = 1.0e-7
for i in range(norbs):
  for j in range(norbs):
      for k in range(norbs):
          for l in range(norbs):
             if(abs(newint2[i,j,k,l]) < 1.0e-14):
                continue;
             if (abs(newint2[i,j,k,l] - newint2[i,l,k,j]) >err 
                 or abs(newint2[i,j,k,l] - newint2[k,j,i,l]) >err
                 or abs(newint2[i,j,k,l] - newint2[k,l,i,j]) >err
                 or abs(newint2[i,j,k,l] - newint2[l,i,j,k]) >err
                 or abs(newint2[i,j,k,l] - newint2[l,k,j,i]) >err
                 or abs(newint2[i,j,k,l] - newint2[j,i,l,k]) >err
                 or abs(newint2[i,j,k,l] - newint2[j,k,l,i]) >err) :
                print(i, j, k, l, newint2[i,j,k,l])
'''



