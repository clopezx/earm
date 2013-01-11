"""
    Overview
    ========
    
    PySB implementations of the apoptosis-necrosis reaction model version 1.0
    (ANRM 1.0) originally published in [Irvin,NECRO2013]_.
    
    This file contains functions that implement the extrinsic apoptosis pathway
    in five modules:
    
    - CD95 Ligation to formation of secondary complex
    - TNFR1 ligation to formation of complex II
    - Secondary complexes to formation of Riptosomes and Necrosomes and Bid
    activation.
    - Execution of apoptosis (MOMP from Albeck)
    - Execution of necrosis
    
    For the (MOMP) segment there are five variants, which correspond to the five
    models described in Figure 11 of [Albeck2008]_:
    
    - "Minimal Model" (Figure 11b, :py:func:`albeck_11b`)
    - "Model B + Bax multimerization" (Figure 11c, :py:func:`albeck_11c`)
    - "Model C + mitochondrial transport" (Figure 11d, :py:func:`albeck_11d`)
    - "Current model" (Figure 11e, :py:func:`albeck_11e`)
    - "Current model + cooperativity" (Figure 11f, :py:func:`albeck_11f`)
    """

import numpy
from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components
from shared_anrm import *
#from earm.shared import *

Model()

def CD95_to_SecondaryComplex_monomers():
    """ Declares Fas ligand, CD95, FADD, Flip_L, Flip_S procaspase8 and Caspase 8.
        
    'bf' is the site to be used for all binding reactions.
        
    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """
    Monomer('Fas', ['blig'])    #Fas Ligand
    Monomer('CD95', ['blig', 'bfad'])           #Fas Receptor (CD95)
    Monomer('FADD', ['brec', 'bc81','bc82'])    #FADD
    Monomer('flip_L', ['bfad'])   #c-Flip[L] binds FADD at bca1 or bca2
    Monomer('flip_S', ['bfad'])   #c-Flip[S] binds FADD at bca1 or bca2
    Monomer('proC8', ['bfad'])            #procaspase 8 binds FADD at bca1 or bca2
    Monomer('C8', ['bf'])                       #active caspase 8

def CD95_to_SecondaryComplex():
    """Defines the interactoins from CD95 ligation of generation of secondary
    complexes as per ANRM 1.0.
    
    Uses Fas, CD95, FADD, flip_L, flip_S, proC8 monomers and C8 active dimers and
    there assicated parameters to generate rules that describe the ligand/receptor
    binding, FADD recruitment, proC8 and c-flip recruitment, activation of caspase
    and release of the secondary complex. This model assumes that one copy of proC8
    binds FADD before c-flip. 
    
    This model converts proC8:proC8 to C8 (active caspase 8 dimer)
    This model also produces Secondary complex, FADD:proC8:c-Flip.
    """
    
    Parameter('KF', 1e-6)
    Parameter('KR', 1e-3)
    Parameter('KC', 1)
    
    Parameter('KF1', 0)
    Parameter('KR1', 1e-1)
    Parameter('KC1', 1)

    Parameter('Fas_0'   ,   3000) # 3000 correspons to 50ng/ml Fas(?)
    Parameter('CD95_0'  ,    200) # 200 receptors per cell
    Parameter('FADD_0'  ,  1.0e3) # molecules per cell (arbitrarily assigned)
    Parameter('flip_L_0',  1.0e2) # molecules per cell
    Parameter('flip_S_0',  1.0e2) # molecules per cell
    Parameter('proC8_0' ,  2.0e4) # procaspase 8 molecules per cell
    Parameter('C8_0'    ,      0) # active caspase 8 dimers per cell.

    Initial(Fas(blig=None), Fas_0)       #Fas Ligand
    Initial(CD95(blig=None, bfad=None), CD95_0)     #Fas Receptor (CD95)
    Initial(FADD(brec=None, bc81=None, bc82=None), FADD_0) #FADD
    Initial(flip_L(bfad=None), flip_L_0)   #c-Flip[L]
    Initial(flip_S(bfad=None), flip_S_0)   #c-Flip[S]
    Initial(proC8(bfad=None), proC8_0)    #procaspase 8
    Initial(C8(bf=None), C8_0)       #caspase 8
   
    # =========================================
    # CD95 ligation and formation of DISC rules
    # -----------------------------------------
    #   Fas + CD95 <-> Fas:CD95
    #   Fas:CD95 + FADD <-> Fas:CD95:FADD
    #   Fas:CD95:FADD + proC8 <-> Fas:CD95:FADD:proC8
    
    #   Fas:CD95:FADD:proC8 + proC8 <-> Fas:CD95:FADD:proC8:proC8 -> Fas:CD95:FADD + C8
    #   Fas:CD95:FADD:proC8 + flip_L <-> Fas:CD95:FADD:proC8:flip_L
    #   Fas:CD95:FADD:proC8 + flip_S <-> Fas:CD95:FADD:proC8:flip_S
    
    #   Fas:CD95:FADD:proC8:flip_L <-> Fas:CD95 + FADD:proC8:flip_L
    #   Fas:CD95:FADD:proC8:flip_S <-> Fas:CD95 + FADD:proC8:flip_S
    # ------------------------------------------

    # -------------DISC assembly----------------
    bind(Fas(blig=None), 'blig',  CD95(blig = None, bfad = None), 'blig', [KF, KR])
    bind(CD95(blig = ANY, bfad = None), 'bfad', FADD(brec = None, bc81 =None, bc82 = None), 'brec', [KF, KR])
    
    bind(FADD(brec = ANY, bc81 = None, bc82 = None),'bc81', proC8(bfad = None), 'bfad', [KF, KR])
    # bind(FADD(brec= None, bc81 = None, bc82 = None), 'bc81', proC8(bfad = None), 'bfad', [KF1, KR1])
    # For simplicity allow proC8 to bind FADD before any c-Flip do.
    bind(FADD(brec = ANY, bc82 = None, bc81 = ANY), 'bc82', flip_L(bfad = None), 'bfad', [KF, KR])
    bind(FADD(brec = ANY, bc82 = None, bc81 = ANY), 'bc82', flip_S(bfad = None), 'bfad', [KF, KR])
    
    
    # procaspase 8 dimerization and activation
    bind(FADD(brec = ANY, bc82 = None, bc81 = ANY), 'bc82', proC8(bfad = None), 'bfad', [KF, KR])
    DISC_proC8 = CD95(blig=ANY, bfad=ANY) % Fas(blig=ANY) % FADD(brec=ANY, bc81=ANY, bc82=ANY) % proC8(bfad=ANY)%proC8(bfad=ANY)
    DISC = CD95(blig=ANY, bfad=ANY) % Fas(blig=ANY) % FADD(brec=ANY, bc81=ANY, bc82=ANY)
    Rule('C8_activation', DISC_proC8 >> DISC + C8(bf=None), KC)
    SecComp = FADD(brec=None, bc81=ANY, bc82=ANY) % proC8(bfad=ANY)%proC8(bfad=ANY)
    Rule('C8_activation2', SecComp >> FADD(brec=None, bc81=None, bc82=None) + C8(bf=None), KC)

    # release of secondary complex from the DISC
    bind(FADD(brec = None, bc82 = ANY, bc81 = ANY), 'brec', CD95(bfad=None), 'bfad', [KF1,KR])

    
#def TNFR1_toSecondaryComplex_monomers()

CD95_to_SecondaryComplex_monomers()
CD95_to_SecondaryComplex()

Observable('ObsFas',  Fas(blig=None))
Observable('ObsFas_CD95', CD95(blig=ANY))
Observable('ObsFas_CD95_FADD', Fas(blig=ANY)%CD95(blig=ANY,bfad=ANY)%FADD(brec=ANY))
Observable('ObsCD95_FADD', CD95(blig=None,bfad=ANY)%FADD(brec=ANY))
Observable('ObsCD95_FADD_proC8', CD95(blig=None, bfad=ANY) % FADD(brec=ANY, bc81=ANY, bc82=None) % proC8(bfad=ANY))
Observable('ObsproC8', proC8(bfad=None))
Observable('ObsC8',  C8(bf=None))
Observable('ObsFADD_flipS_proC8', FADD(brec=None, bc81=1, bc82=2) % flip_S(bfad=2) % proC8(bfad=1))

#from pysb.integrate import odesolve
#m = Model
#t = numpy.linspace(0,30000,100)
#y = odesolve(m,t)

#plot(t,ObsFas)
