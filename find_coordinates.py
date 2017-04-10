#!/usr/bin/env python

# Autor : Alexandre SImoneau
# Date : 2015/10/25
# Code pour trouver les points les plus proches de ceux choisis pour les zones
# en utilisant le fichier de conversion fourni avec les intrants.
#
# Il faut produire un fichier contenant les parametres des zones.
# Les lignes qui commencent avec '#' ne sont pas lues.
# Ex :
#
# # longitude   latitude    rayon
# -75.12345   45.12345    15
# ...
#
# Produit un fichier appele 'zones.dat' contenant les parametres des zones.
# Ces valeurs de X et Y vont etre a corriger pour correspondre aux valeurs
# attendues par ILLUMINA. La formule de correction est disponible sur le wiki.

# importation de numpy : necessaire pour pouvoir utiliser ses fonctions
import numpy as np

# importer les donnees de conversions fourines
data = np.loadtxt("row_column_to_x_y_lon_lat.csv",skiprows=1,delimiter=",")
# importer les parametres des zones trouvees avec le site web
# fichier au format suivant :
# lon   lat   rayon
zones = np.loadtxt("zones_coords.dat",ndmin=2)

# extraction des coordonnes
coords = zones[:,:2]
lonlat = data[:,-2:]

# calcul du carre de la distance entre chaque coordonnee donnee et chaque coordonnee dans la grille
# min(x^2) = min(x)
dist_2 = np.sum((lonlat - coords[:,np.newaxis])**2, axis=-1)

# on trouve les positions des points les plus proches dans la liste
index = np.argmin(dist_2, axis=-1)
# on trouve les coordonnee des points trouves
pos = data[index,:2]

# on ajoute les rayons des zones a leurs coordonnes
new_zones = np.concatenate( [ pos, zones[:,-1,np.newaxis] ], axis=-1 )

# on enregistre dans un fichier
np.savetxt("zones.dat", new_zones, "%g")

# On corrige pour ILLUMINA
H = int(raw_input("Quelle est la hauteur de l'image ? (px) : "))

new_zones[:,0] += 1 # X -> X+1
new_zones[:,1] = H-new_zones[:,1] # Y -> H-Y

np.savetxt("zones_corr.dat", new_zones, "%g")



