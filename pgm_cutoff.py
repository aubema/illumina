cutoff = 256

with open("pgms/stable_lights.pgm") as file:
	data = file.read().split('\n')[:-1]

header = data[:8]
data = map(lambda s: '0' if int(s)<cutoff else s, data[8:])

with open("pgms/stable_lights_lowcut.pgm",'w') as file:
	file.write('\n'.join(header+data))
