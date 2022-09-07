# This code takes a plot file as input, averages the variables chosen in 'varPos' in z direction,
# generates a plot file with the averaged values, and makes a list of averaged values. Each line
# of this new list contains: x, y, averaged value of variables appeared in 'varPos'.

from numpy import *

varPos = [4,5,6,7,8,9,10,11,12]       # <<<--------- position of variables in the plot file:

varPos = [x-1 for x in varPos]

# position of ariables in the list, data:
variables = range(len(varPos))

data = []
for i in variables: data.append({})

moviename = 'movie0002'

fold = open(moviename+'.plt','r')

# read in and construct numerator sums
print 'reading data...'
for line in fold:
   if ('TITLE' not in line) and ('TEXT' not in line) and ('VARIABLES' not in line) and ('ZONE' not in line):
      (x,y) = (float(line.split()[0]),float(line.split()[1]))
      
      for var in variables:
         if (x,y) not in data[var]:
            data[var][(x,y)] = float(line.split()[varPos[var]])
            data[var][(x,y,1)] = 1
         else:
            data[var][(x,y)] = data[var][(x,y)]+float(line.split()[varPos[var]])
            data[var][(x,y,1)] = data[var][(x,y,1)] + 1
fold.close()

# construct averages based on npz
print 'constructing averages...'
for var in variables:
   for key in data[var]:
      keya = key
      keyb = (key[0],key[1],1)
      if keya != keyb:
         data[var][keya] = data[var][keya]/data[var][keyb]

fold = open(moviename+'.plt','r')
fnew = open(moviename+'.tec','wb')
output = open('list.txt','w')

# write new file with averaged data
print 'writing data...'
for line in fold:
   if ('TITLE' in line) or ('TEXT' in line) or ('VARIABLES' in line) or ('ZONE' in line):
      fnew.write(line)
   else:
      (x,y) = (float(line.split()[0]),float(line.split()[1]))
      z = float(line.split()[2])
      
      if (x,y) in data[0]:
         newline = line.split()
         for var in variables: newline[varPos[var]] = str(data[var][(x,y)])
         for item in newline:
            fnew.write(str(item)+' ')
         fnew.write('\n')
         
         if z<1.0e-3:
            output.write(str(x).ljust(18)+', ')
            output.write(str(y).ljust(18)+', ')
            for var in variables[:-1]:
               output.write(str(data[var][(x,y)]).ljust(18)+', ')
            var = variables[-1]
            output.write(str(data[var][(x,y)]).ljust(18))
            output.write('\n')
         
      else:
         print 'we have a missing value from the dictionary!'
         
fnew.close()
fold.close()
output.close()
