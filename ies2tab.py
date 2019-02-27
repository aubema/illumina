#!/usr/bin/python2
# cedi est un commentaire
import numpy as np

filename = raw_input("IES filename : ")
with open(filename) as f:
	data = f.readlines()
data = filter(None,map(str.strip, data))
data = filter(lambda s: s[0].isdigit(), data)
data = sum(map(lambda s: s.split(),data),[])

nbV  = int(data[3])
nbH  = int(data[4])
Vang = data[13:13+nbV]
Hang = data[13+nbV:13+nbV+nbH]
data = map(None,*[iter(data[13+nbV+nbH:])]*nbH)

outname = raw_input("out filename : ")
with open(outname,'w') as f:
	f.write("%d %d\n" % (nbV,nbH))
	f.write(' '.join(Vang)+'\n')
	f.write(' '.join(Hang)+'\n')
	f.write('\n'.join([' '.join(d) for d in data])+'\n')
