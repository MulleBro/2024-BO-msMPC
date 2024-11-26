# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 11:13:27 2023

@author: TobiasBrockhoff
"""

import time
import do_mpc
import numpy as np
import pdb
import sys
import os
rel_do_mpc_path = os.path.join('..', '..')
sys.path.append(rel_do_mpc_path)


def model_names(model):
    x_names=[]
    for i in range(model.n_x):
        xi=model._x.getLabel(i)
        xi=xi.replace('[','')
        xi=xi.replace(',0]','')
        x_names.append(xi)
        
    aux_names=[]
    for i in range(model.n_aux):
        auxi=model._aux.getLabel(i)
        auxi=auxi.replace('[','')
        auxi=auxi.replace(',0]','')
        aux_names.append(auxi)
    p_names=[]
    for i in range(model.n_p):
        pi=model._p.getLabel(i)
        pi=pi.replace('[','')
        pi=pi.replace(',0]','')
        p_names.append(pi)
        
    u_names=[]
    for i in range(model.n_u):
        ui=model._u.getLabel(i)
        ui=ui.replace('[','')
        ui=ui.replace(',0]','')
        u_names.append(ui)

    return x_names,p_names, aux_names, u_names
        
    
def model_list(model):
    print('STATES _x: ')
    for i in range(model.n_x):
        print(str(i)+' = '+model._x.getLabel(i))

    print('INPUTS _u: ')
    for j in range(model.n_u):
        print(str(j)+' = '+model._u.getLabel(j))

    print('AUX _aux: ')
    for k in range(model.n_aux):
        print(str(k)+' = '+model._aux.getLabel(k))

    print('PARAMETERS _p: ')
    for l in range(model.n_p):
        print(str(l)+' = '+model._p.getLabel(l))

    print('PARAMETERS _tvp: ')
    for l in range(model.n_tvp):
        print(str(l)+' = '+model._tvp.getLabel(l))

#getting index of state -> more robust different models
def get_ind_x(model,state):
    search_str='['+state+',0]'
    found=False
    
    for i in range(model.n_x):
        if model._x.getLabel(i)==search_str:
            found=True
            break
    if found==False:
        i='STATE NOT FOUND IN MODEL'
    return i

def get_ind_p(model,param):
    search_str='['+param+',0]'
    found=False
    
    for i in range(model.n_p):
        if model._p.getLabel(i)==search_str:
            found=True
            break
    if found==False:
        i='PARAM NOT FOUND IN MODEL'
    return i

def get_ind_aux(model,aux):
    search_str='['+aux+',0]'
    found=False
    
    for i in range(model.n_aux):
        if model._aux.getLabel(i)==search_str:
            found=True
            break
    if found==False:
        i='AUX NOT FOUND IN MODEL'
    return i



