import struct
import numpy as np
import scipy as sp
import matplotlib.pyplot as m


#simple single run example, add a for loop for parameter looping, etc. 




fileName="./data/testdata.dat_5"
with open(fileName, mode='rb') as file:
	dataBuffer=file.read(8)
	steps=struct.unpack('d',dataBuffer)
	steps=int(steps[0]);	
	dataBuffer=file.read(8)
	neutrons=struct.unpack('d',dataBuffer)
	neutrons=int(neutrons[0]);
	dataBuffer=file.read(8)
	variables=struct.unpack('d',dataBuffer)
	variables=int(variables[0]);
	dataBuffer=file.read()
	#the next line turns out to be pretty slow in python, 
	#I'm thinking about adding option of averageing in simulation. 
	data=struct.unpack("d"*(steps*neutrons*variables),dataBuffer)

print 'There are ' + str(neutrons) + ' neutrons for ' +str(steps) + ' time steps, each with ' + str(variables) + ' variables'

data=np.asanyarray(data)
#remember python is [start:stop:step] (stupid)

data=np.reshape(data,(neutrons,steps,variables))


sx=data[0:neutrons,0:steps-1,0]
sy=data[0:neutrons-1,0:steps-1,1]
sz=data[0:neutrons-1,0:steps-1,2]
phase=data[0:neutrons-1,0:steps-1,3]
xpos=data[0:neutrons-1,0:steps-1,4]
ypos=data[0:neutrons-1,0:steps-1,5]
zpos=data[0:neutrons-1,0:steps-1,6]
xvel=data[0:neutrons-1,0:steps-1,7]
yvel=data[0:neutrons-1,0:steps-1,8]
zvel=data[0:neutrons-1,0:steps-1,9]
t=data[0,0:steps-1,10]

if neutrons>1:
	Sx=np.mean(sx,axis=0)
	Sy=np.mean(sy,axis=0)
	Sz=np.mean(sz,axis=0)
	Phase=np.mean(phase,axis=0)
else:
	Sx=sx;
	Sy=sy;
	Sz=sz;
	Phase=phase;
#print str(Phase)

fig, ax = m.subplots()

#single neutron, 
#ax.plot(t,phase[0,0:steps-1])
#z position plot, 
#ax.plot(t,xpos[0,0:steps-1])
#ax.plot(t,ypos[0,0:steps-1])
#ax.plot(t,zpos[0,0:steps-1])
#mean phase plot. 
ax.plot(t,Sx)

#mean Sx plot. 
#ax.plot(t,Sx)


ax.set(xlabel='time (s)', ylabel='phase',
        title='neutron precession')
ax.grid()

fig.savefig("Sphere_mirror_shift_x_spin.png")
m.show()