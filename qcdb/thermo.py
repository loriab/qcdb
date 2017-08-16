#
# @BEGIN LICENSE
#
# QCDB: quantum chemistry common driver and databases
#
# Copyright (c) 2011-2017 The QCDB Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of QCDB.
#
# QCDB is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# QCDB is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with QCDB; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import numpy as np
from numpy import linalg as la
import math
import qcdb
from periodictable import *
import psi4

pi = np.pi

def mass_weight_Hessian(Hessian, masses):
    # Hessian is numpy matrix, and masses is numpy array
    masses_mat = np.zeros(3*len(masses))
    for i in range(len(masses)):
        masses_mat[i*3:i*3+3] = masses[i]
   
    masses_mat = 1/(np.sqrt(masses_mat))
    masses_mat = masses_mat*np.identity(3*len(masses)) 

    mwHessian = np.dot(np.dot(masses_mat,Hessian),masses_mat)
    
    return mwHessian   

def get_rot_const(molecule, text):
    #molecule is qcdb molecule object
    #get rotational constants from qcdb
    rot_const = np.array(molecule.rotational_constants())
    
    # moments in amu * cm^2
    moments = []
    for i in rot_const:
        if i is None: moments.append(i)
        else: moments.append(qcdb.physconst.psi_h/(i * 8. * pi * pi * qcdb.physconst.psi_c))

    text+="\n ==> Moments of Inertia <==\n"
    text+="%10s  %10s\n" % ('','amu * cm^2')
    if moments[0] is None: text+="%10s  %10s\n" % ('A:','*')
    else: text+="%10s  %10.5e\n" % ('A:',moments[0])
    if moments[1] is None: text+="%10s  %10s\n" % ('B:','*')
    else: text+="%10s  %10.5e\n" % ('B:',moments[1])
    if moments[2] is None: text+="%10s  %10s\n" % ('C:','*')
    else: text+="%10s  %10.5e\n" % ('C:',moments[2])
 
    rot_conv = (100.*qcdb.physconst.psi_c)/1.e9

    text+="\n ==>Rotational Constants <==\n"
    text+="%10s  %10s  %10s\n" % ('','cm^-1','GHz')
    if rot_const[0] is None: text+="%10s  %10s  %10s\n" % ('A:','*','*')
    else:text+="%10s  %10.5f  %10.5f\n" % ('A:',rot_const[0],rot_conv*rot_const[0])
    if rot_const[1] is None:text+="%10s  %10s  %10s\n" % ('B:','*','*')
    else:text+="%10s  %10.5f  %10.5f\n" % ('B:',rot_const[1],rot_conv*rot_const[1])
    if rot_const[2] is None:text+="%10s  %10s  %10s\n" % ('C:','*','*')
    else:text+="%10s  %10.5f  %10.5f\n" % ('C:',rot_const[2],rot_conv*rot_const[2])
   
    return rot_const, text

def get_point_group_info(molecule, text):
    #molecule is qcdb molecule object
    # grab pg
    pg = str(molecule.full_pg).split(':')[-1].upper()
    # grab full pg number
    pg_n = molecule.full_pg_n()
    # rotor type 
    rot_type = molecule.rotor_type()

    rot_symm_type = {"ATOM":1,"C1":1,"CI":1,"CS":1,"C_INF_V":1,"D_INF_H":2,"TD":12,"OH":24,"IH":60,
    "CN":pg_n,"CNV":pg_n,"CNH":pg_n,"DN":2.*pg_n,"DND":2.*pg_n,"DNH":2.*pg_n,"SN": pg_n / 2.}
    
    if pg in rot_symm_type: rot_symm_num = rot_symm_type[pg]
    else: raise qcdb.exceptions.ValidationError("thermo_analysis(): Could not interpret molecular point group.")    

    return pg, pg_n, rot_symm_num, rot_type, text
    
def get_vibfreq_and_temp(w,rot_type, text): 
    #w is numpy array of vibrational modes 
    text+='\n ==>  Harmonic Vibrational Frequencies (cm^-1)  <==\n'
    k_convfactor = qcdb.physconst.psi_hartree2J/(qcdb.physconst.psi_bohr2m * qcdb.physconst.psi_bohr2m * qcdb.physconst.psi_amu2kg)
    cm_convfactor = 1.0/(2.0 * pi * qcdb.physconst.psi_c * 100.) 
    
    if rot_type == 'RT_ATOM': minus = len(w)
    elif rot_type == 'RT_LINEAR': minus = 5
    else: minus = 6

    #len(w) - minus is number of vibrational freqs
    sw = sorted(w, key=abs, reverse = True)
    modes = sorted(sw[:(len(w)-minus)])
    nonmodes = sorted(sw[(len(w)-minus):])
    #rotational and translational freq removed from vibfreq
    vibfreq = []
    m_remove = []
    for i in modes: 
        if abs(i) < 3.78435124825e-14:
            nonmodes.append(i)
            nonmodes = sorted(nonmodes)
            m_remove.append(i)
            if rot_type == 'RT_LINEAR': 
                text+="\nWARNING: Degenerate bending mode in near linear molecule\n"
                text+="not found. Treating thermodynamic calculations as ASYMMETRIC_TOP.\n\n"
                rot_type = 'RT_ASYMMETRIC_TOP'
                minus = 6
        elif i < 0.0: 
            vibfreq.append(-1.*cm_convfactor * pow(abs(k_convfactor * i),0.5))
        else: 
            vibfreq.append(cm_convfactor * pow(abs(k_convfactor * i), 0.5))
    for i in m_remove: modes.remove(i)
    garbagefreq = []
    for i in nonmodes:
        if i < 0.0: garbagefreq.append(-1.*cm_convfactor * pow(abs(k_convfactor * i),0.5))
        else: garbagefreq.append(cm_convfactor * pow(abs(k_convfactor * i), 0.5))

    #rotational and translational freq removed from vibtemp
    vibtemp = []
    for i in range(len(modes)):
        vibtemp.append(100 * qcdb.physconst.psi_h * qcdb.physconst.psi_c * vibfreq[i] / qcdb.physconst.psi_kb)

    garbagetemp = []
    for i in range(len(nonmodes)): 
        garbagetemp.append(100 * qcdb.physconst.psi_h * qcdb.physconst.psi_c * garbagefreq[i] / qcdb.physconst.psi_kb)
    
    text+= "%20s  %13s\n" % ('Vib. Freq. (cm^-1)', 'Vib. Temp (K)')
    for i in range(len(vibfreq)):
        text+= "%20.3f  %13.3f\n" % (vibfreq[i],vibtemp[i])
    for i in range(len(garbagefreq)):
        text+= "%20.3f  %13.3f  %30s\n" % (garbagefreq[i],garbagetemp[i],'[Translational or Rotational]')
    
    #return rot_type in case it changes for near linear molecule
    return vibfreq, vibtemp, rot_type, text

def compute_therm_terms(E0, molecule, vibfreq, vibtemp, rot_symm_num, rot_const, rot_type, T, P, mw, text):
    #E0 is current energy (float), molecule is qcdb molecule
    #object, vibfreq is numpy array, vibtemp is numpy array,
    #rot_symm_num is integer, rot_const is numpy array,
    #rot_type is string, T is float, P is float, mw is float 

    terms = {'Electronic':0,'Translational':1,'Rotational':2,'Vibrational':3,'Total':4}
    kterms = ['Electronic','Translational','Rotational','Vibrational','Total']

    E = {} 
    Cv = {} 
    Cp = {} 
    S = {}
    H = {}
    G = {}
     
    beta = 1./(qcdb.physconst.psi_kb*T)

    if rot_type == "RT_ATOM": 
        E['Rotational'] = Cv['Rotational'] = Cp['Rotational'] = S['Rotational'] = 0.0
    elif rot_type == "RT_LINEAR": 
        E['Rotational'] = T
        Cv['Rotational'] = Cp['Rotational'] = 1.0
        q_rot = 1. / (beta * rot_symm_num*100*qcdb.physconst.psi_c*qcdb.physconst.psi_h*rot_const[1])
        S['Rotational'] = 1.0 + math.log(q_rot) 
    else: # EVERYTHING ELSE
        E['Rotational'] = 1.5*T
        Cv['Rotational'] = Cp['Rotational'] = 1.5
        phi_A, phi_B, phi_C = rot_const*100*qcdb.physconst.psi_c*qcdb.physconst.psi_h/qcdb.physconst.psi_kb
        q_rot = math.sqrt(pi)*pow(T,1.5)/(float(rot_symm_num)*math.sqrt(phi_A*phi_B*phi_C))
        S['Rotational'] = 1.5 + math.log(q_rot)
    H['Rotational'] = E['Rotational'] 
    G['Rotational'] = H['Rotational'] - T * S['Rotational']
    
    ZPVE = 0.0
    E['Vibrational'] = Cv['Vibrational'] = Cp['Vibrational'] = S['Vibrational'] = 0.0
    for i in range(len(vibfreq)):
        if vibtemp[i] < 0.0: 
            text+= "WARNING: vibration with imaginary frequency neglected in vibrational contributions.\n"
        else:
            if vibtemp[i] < 900.0: text+="WARNING: used thermodynamic relations are not appropriate for low frequency modes.\n"
            rT = vibtemp[i]/T
            E['Vibrational'] += vibtemp[i]*(0.5+1.0/(math.exp(rT)-1.))
            Cv['Vibrational'] += math.exp(rT)*math.pow((rT/(math.exp(rT)-1.)),2.)
            Cp['Vibrational'] += math.exp(rT)*math.pow((rT/(math.exp(rT)-1.)),2.)
            S['Vibrational'] +=  rT/(math.exp(rT)-1.) - math.log(1.-math.exp(-rT))
            ZPVE += vibfreq[i]/2.
            H['Vibrational'] = E['Vibrational']
            G['Vibrational'] = H['Vibrational'] - T * S['Vibrational']
    # Electronic
    q_elec = molecule.multiplicity()
    H['Electronic'] = 0.0
    E['Electronic'] = 0.0
    Cp['Electronic'] = 0.0
    Cv['Electronic'] = 0.0
    S['Electronic'] = math.log(q_elec)
    G['Electronic'] = H['Electronic'] - T * S['Electronic']
    # Translational
    E['Translational'] = 1.5 * T
    Cv['Translational'] = 1.5
    Cp['Translational'] = 2.5
    q_trans = pow(2.0*pi*mw*qcdb.physconst.psi_amu2kg/(beta*qcdb.physconst.psi_h*qcdb.physconst.psi_h),1.5)*qcdb.physconst.psi_na/(beta*P)
    S['Translational'] = (5.0/2.0) + math.log(q_trans/qcdb.physconst.psi_na)
    H['Translational'] = 2.5 * T
    G['Translational'] = H['Translational'] - T * S['Translational']
    
    ## Unit conversions
    R_to_cal = qcdb.physconst.psi_R/qcdb.physconst.psi_cal2J
    for i in E: E[i] *= R_to_cal/1000.0
    for i in H: H[i] *= R_to_cal/1000.0
    for i in G: G[i] *= R_to_cal/1000.0
    for i in Cv: Cv[i] *= R_to_cal
    for i in Cp: Cp[i] *= R_to_cal
    for i in S: S[i] *= R_to_cal
    # Total energy 
    E['Total'] = E['Translational'] + E['Rotational'] + E['Vibrational'] + E['Electronic']
    Cp['Total'] = Cp['Translational'] + Cp['Rotational'] + Cp['Vibrational'] + Cp['Electronic']
    Cv['Total'] = Cv['Translational'] + Cv['Rotational'] + Cv['Vibrational'] + Cv['Electronic']
    S['Total'] = S['Translational'] + S['Rotational'] + S['Vibrational'] + S['Electronic']
    H['Total'] = H['Translational'] + H['Rotational'] + H['Vibrational'] + H['Electronic']
    G['Total'] = G['Translational'] + G['Rotational'] + G['Vibrational'] + G['Electronic']
    
    ZPVE_au = ZPVE * 100 * qcdb.physconst.psi_h * qcdb.physconst.psi_c / qcdb.physconst.psi_hartree2J
    DU = E['Total'] * 1000.0 * qcdb.physconst.psi_cal2J / qcdb.physconst.psi_na / qcdb.physconst.psi_hartree2J
    DH = DU + qcdb.physconst.psi_kb * T / qcdb.physconst.psi_hartree2J
    DG = DH - T * S['Total'] * qcdb.physconst.psi_cal2J / qcdb.physconst.psi_na / qcdb.physconst.psi_hartree2J
   

    text+="\n==> Components <==\n\n"
    text+="Entropy, S\n"
    for i in range(5):
        term = kterms[i]
        ScalmolK = S[term]
        SjmolK = ScalmolK * qcdb.physconst.psi_cal2J
        SmEhK = ScalmolK / qcdb.physconst.psi_hartree2kcalmol
        if kterms[i] == "Electronic": text+="  %15s S  %11.3lf [cal/(mol K)]  %11.3lf [J/(mol K)]  %15.8lf [mEh/K] (multiplicity = %d)\n" % (term,ScalmolK,SjmolK,SmEhK,molecule.multiplicity())
        elif kterms[i] == "Translational": text+="  %15s S  %11.3lf [cal/(mol K)]  %11.3lf [J/(mol K)]  %15.8lf [mEh/K] (mol. weight = %.4f [u], P = %.2f [Pa])\n" % (term,ScalmolK,SjmolK,SmEhK,mw,P)
        elif kterms[i] == "Rotational": text+="  %15s S  %11.3lf [cal/(mol K)]  %11.3lf [J/(mol K)]  %15.8lf [mEh/K] (symmetry no. = %d)\n" % (term,ScalmolK,SjmolK,SmEhK, rot_symm_num)
        elif kterms[i] == "Total": text+=" Total S %10s  %11.3lf [cal/(mol K)]  %11.3lf [J/(mol K)]  %15.8lf [mEh/K]\n" % ('',ScalmolK,SjmolK,SmEhK)
        else: text+="  %15s S  %11.3lf [cal/(mol K)]  %11.3lf [J/(mol K)]  %15.8lf [mEh/K]\n" % (term,ScalmolK,SjmolK,SmEhK)
    
    text+="\nConstant volume heat capacity, Cv\n"
    for i in range(5):
        term = kterms[i]
        CvcalmolK = Cv[term]
        CvjmolK = CvcalmolK * qcdb.physconst.psi_cal2J
        CvmEhK = CvcalmolK / qcdb.physconst.psi_hartree2kcalmol
        if term == "Total": text+=" Total Cv %10s  %11.3lf [cal/(mol K)]  %11.3lf [J/(mol K)]  %15.8lf [mEh/K]\n" % ('',CvcalmolK,CvjmolK,CvmEhK)
        else: text+="  %15s Cv  %11.3lf [cal/(mol K)]  %11.3lf [J/(mol K)]  %15.8lf [mEh/K]\n" % (term,CvcalmolK,CvjmolK,CvmEhK)
    
    text+="\nConstant pressure heat capacity, Cp\n"
    for i in range(5):
        term = kterms[i]
        CpcalmolK = Cp[term]
        CpjmolK = CpcalmolK * qcdb.physconst.psi_cal2J
        CpmEhK = CpcalmolK / qcdb.physconst.psi_hartree2kcalmol
        if term == "Total": text+=" Total Cp %10s  %11.3lf [cal/(mol K)]  %11.3lf [J/(mol K)]  %15.8lf [mEh/K]\n" % ('',CpcalmolK,CpjmolK,CpmEhK)
        else: text+="  %15s Cp  %11.3lf [cal/(mol K)]  %11.3lf [J/(mol K)]  %15.8lf [mEh/K]\n" % (term,CpcalmolK,CpjmolK,CpmEhK)
        

    text+="\n==> Energy Analysis <==\n"
    text+="\nRaw electronic energy, E0\n"
    text+="Total E0, Electronic energy at well bottom at 0 [K] %15s %15.8lf [Eh]\n\n" % ('',E0) 
    text+="Zero-point energy, ZPE_vib = Sum_i nu_i / 2\n"
    for i in range(5):
        ZPVE_kcalmol = ZPVE_au*qcdb.psi_hartree2kcalmol
        ZPVE_kjmol = ZPVE_kcalmol*qcdb.psi_cal2J
        ZPVE_cm = ZPVE_au*qcdb.psi_hartree2wavenumbers
        if kterms[i] == "Total": text+=" Correction ZPE %5s %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]  %15.3lf [cm^-1]\n" % ('',ZPVE_kcalmol,ZPVE_kjmol,ZPVE_au,ZPVE_cm)
        elif kterms[i] == "Vibrational": text+="  %15s ZPE %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]  %15.3lf [cm^-1]\n" % (kterms[i],ZPVE_kcalmol,ZPVE_kjmol,ZPVE_au,ZPVE_cm)
        else: text+="  %15s ZPE %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % (kterms[i],0.0,0.0,0.0)
    text+="Total ZPE, Electronic energy at 0 [K] %29s %15.8lf [Eh]\n\n" % ('',E0+ZPVE_au) 
    
    text+="Thermal Energy, E (includes ZPE)\n"
    for i in range(5):
        term = kterms[i]
        Ekcalmol = E[term]
        Ekjmol = Ekcalmol*qcdb.psi_cal2J
        Eh = Ekcalmol/qcdb.psi_hartree2kcalmol
        if term == "Total": text+=" Correction E %7s %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % ('',Ekcalmol,Ekjmol,Eh)
        elif term == "Vibrational": text+="  %15s E   %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % (term,Ekcalmol,Ekjmol,Eh)
        else: text+="  %15s E   %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % (term,Ekcalmol,Ekjmol,Eh)
    text+="Total E, Electronic energy at %.2f [K] %26s %15.8lf [Eh]\n\n" % (T,'',E0+DU) 
    
    text+="Enthalpy, H_trans = E_trans + k_B * T\n"
    for i in range(5):
        term = kterms[i]
        Hkcalmol = H[term]
        Hkjmol = Hkcalmol*qcdb.psi_cal2J
        Hh = Hkcalmol/qcdb.psi_hartree2kcalmol
        if term == "Total": text+=" Correction H %7s %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % ('',Hkcalmol,Hkjmol,Hh)
        elif term == "Vibrational": text+="  %15s H   %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % (term,Hkcalmol,Hkjmol,Hh)
        else: text+="  %15s H   %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % (term,Hkcalmol,Hkjmol,Hh)
    text+="Total H, Enthalpy at %.2f [K] %35s %15.8lf [Eh]\n\n" % (T,'',E0+DH) 

    text+="Gibbs free energy, G = H - T * S\n"
    for i in range(5):
        term = kterms[i]
        Gkcalmol = G[term]
        Gkjmol = Gkcalmol*qcdb.psi_cal2J
        Gh = Gkcalmol/qcdb.psi_hartree2kcalmol
        if term == "Total": text+=" Correction G %7s %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % ('',Gkcalmol,Gkjmol,Gh)
        elif term == "Vibrational": text+="  %15s G   %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % (term,Gkcalmol,Gkjmol,Gh)
        else: text+="  %15s G   %11.3lf [kcal/mol]  %11.3lf [kJ/mol]  %15.8lf [Eh]\n" % (term,Gkcalmol,Gkjmol,Gh)
    text+="Total G, Free enthalpy at %.2f [K] %30s %15.8lf [Eh]\n\n" % (T,'',E0+DG) 
    
    return text  
 
def do_analysis(wfn, E0, Hessian, molecule, T, P, text):
    #E0 is current energy (float), Hessian is numpy matrix, molecule is qcdb 
    #molecule object, T is float, and P is float.
    #get point group info from qcdb
    pg, pg_n, rot_symm_num, rot_type, text  = get_point_group_info(molecule, text)

    #print parameters
    text+="\nData used to determine thermochemical information:\n\n" 
    text+="%30s %10s\n"%("Temperature (K):",str(T))
    text+="%30s %10s\n"%("Pressure (Pa):",str(P))
    text+="%30s %10s\n"%("Multiplicity:",molecule.multiplicity())
    text+="%30s %10s\n"%("Full Point Group:",pg)  
    text+="%30s %10s\n"%("Rotor Type:",rot_type[3:])
    text+="%30s %10s\n"%("Rotational Symmetry Number:",rot_symm_num)
   
    pg = molecule.find_point_group()
    
    mints = psi4.MintsHelper(wfn.basisset()) 
    cdsalcs = mints.cdsalcs(0xFF, True, True)
    U = np.asarray(cdsalcs.matrix())
   
    k_convfactor = qcdb.physconst.psi_hartree2J/(qcdb.physconst.psi_bohr2m * qcdb.physconst.psi_bohr2m * qcdb.physconst.psi_amu2kg)
    cm_convfactor = 1.0/(2.0 * pi * qcdb.physconst.psi_c * 100.) 
    count = 0
 
    text+= "\n ==> Nuclear Masses <== \n"

    masses = []

    for i in range(molecule.natom()):
        masses.append(molecule.mass(i))

    for i in range(len(masses)):
        text+= "%5s:  %10.5f\n"%(qcdb.periodictable.z2el[molecule.Z(i)],masses[i])
    
    #get rotational constants from qcdb
    rot_const, text = get_rot_const(molecule, text)
    
    mw = sum(masses)
    # get new mass-weighted hessian using current masses
    mwHessian = mass_weight_Hessian(Hessian, masses)
    blocked_H = np.dot(np.dot(U,mwHessian),U.T)
    w, v = np.linalg.eig(blocked_H)

    # get eigenvalues, eigenvectors of mass-weighted hessian 
    w, v = la.eig(mwHessian)

    w = np.real(w)

    # get vibrational freqs and temps, pass rot_type in case it changes
    #for near linear molecules
    vibfreq, vibtemp, rot_type, text = get_vibfreq_and_temp(w, rot_type, text)
    
    #Fail if not all irreps are found
    #if len(vibfreq) != 

    # compute thermodynamic terms
    text = compute_therm_terms(E0, molecule, vibfreq, vibtemp, rot_symm_num, rot_const, rot_type, T, P, mw, text)

    return text

def print_header(text):

    text+='\n\n**********************************************\n'
    text+='*            Thermodynamic Analysis          *\n'
    text+='*      B. Bakr, A. Schile, & L. Burns 2016   *\n'
    text+='**********************************************\n\n'

    return text

def thermo_analysis(wfn, molecule, E0, Hessian, isotopes = [], T = 298.15, P = 101325.00):
    """Conduct thermodynamic analysis from a given molecule, energy,
    and Hessian. Molecular coordinates are assumed to be in Angstrom
    and are passed as a string. Energy is is assumed to be in atomic
    units. Hessian is assumed to be in atomic units / Bohr^2. Isotopic
    masses to be subsituted for standard atomic masses may be passed
    in as list of dictionaries. Temperature is set to 298.15 Kelvin by
    default, and Pressure is set to 101315.00 Pascals by default."""

    text = ''

    #makes molecule of molecule so that masses can be manipulated for isotopic studies
    molecule = qcdb.Molecule(molecule)
    molecule.update_geometry()
 
    #if no isotopic subsittutions given, do standard masses
    if len(isotopes) == 0: isotopes.append({})    

    #for each set of isotopic subs
    for dic in isotopes: 
        text = print_header(text)
        for key in dic.keys():
            newmass = 0.0
            if type(dic[key]) is str: newmass = eliso2mass[dic[key]]
            else: newmass = float(dic[key])
            text+="Changing mass of atom number %s (%s) from %.8f amu to %.8f amu.\n" % (str(key), qcdb.periodictable.z2el[molecule.Z(key-1)], molecule.mass(key-1), newmass)
            #replace mass of atom key-1 with user-defined mass
            molecule.set_mass(key-1, newmass)
        
        text = do_analysis(wfn, E0, Hessian, molecule, T, P, text)
        #reset to standard masses
        for key in dic.keys():
            molecule.set_mass(key-1, el2mass[molecule.symbol(key-1)])
   
    return text
