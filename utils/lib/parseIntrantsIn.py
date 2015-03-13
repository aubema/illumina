# -*- coding: utf-8 -*-
"""Évaluation de la configuration des intrants
"""
from  ConfigParser import ConfigParser
from datetime import datetime
import logging

class Domaine:
    """
    Caractéristique du domaine
    Ces caractéristiques comprennent:
    - le rectangle englobant
    - la taille du pixel
    - la date de modélisation
    - la longueur d'onde
    
    Consulter le fichier intrants.in pour plus d'information.
    """
    def __init__(self,config):
        self._date_str = config.get('domaine','date')
        self.date = datetime.strptime( self._date_str, "%Y-%m-%d")
        self.latMin = float(config.get('domaine', 'latMin'))
        self.latMax = float(config.get('domaine', 'latMax'))
        self.lonMin = float(config.get('domaine', 'lonMin'))
        self.lonMax = float(config.get('domaine', 'lonMax'))
        self._check_bbox()
        
    def _check_bbox(self):
        if self.latMax <= self.latMin:
            raise InputError('Latitude maximale inférieure ou égale à la latitude minimale')
        if self.lonMax <= self.lonMin:
            raise InputError('Longitude maximale inférieure ou égale à la longitude minimale')
            
        

class Observatur:
    "Caractéristique de l'observateur"

class Meteo:
    "Caractéristique météorologique et climatique"

class Intrants:
    "Description des intrants pour un modélisation"
    def __init__(self,file):
        config = ConfigParser()
        read_file = config.read(file)
        if len(read_file) == 0:
            file.seek(0)
            config.readfp(file)
        self.domaine = Domaine(config)
        
class InputError(Exception):
    """Exception lorsqu'une erreur logique est détectée dans les intrants
    
    Attributs:
        expr -- expression où l'erreur est détectée
        msg -- message d'erreur
        
    """
        
    def __init__(self, msg):
        self.msg = msg