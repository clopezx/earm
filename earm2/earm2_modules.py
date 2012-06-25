from pysb import *
from pysb.util import alias_model_components
from earm2.macros import *

def embedded():
    alias_model_components()

    # tBID to MOMP 
    # ======================
    # tBid and Bad spontaneously insert into the mito membrane
    # --------------------------------------------------------
    translocate(Bid(bf=None), 'T', 'M')
    translocate(Bad(bf=None), 'C', 'M')

    # tBid in the mitochondria activates Bax(C) and Bak(M)
    #-----------------------------------------------------
    bid_effector_rates = [1e-7, 1e-3, 1]
    catalyze(Bid(state='M'), Bax(state='C'), Bax(state='A'), bid_effector_rates)
    catalyze(Bid(state='M'), Bak(state='M'), Bak(state='A'), bid_effector_rates)

    # Autoactivation, Bax and Bak activate their own kind, but only when
    # free (i.e. not part of a pore complex)
    # ------------------------------------------------------------------
    effector_auto_rates = [1e-7, 1e-3, 1]
    catalyze(Bax(state='A', s1=None, s2=None), Bax(state='C'), Bax(state='A'), effector_auto_rates)
    catalyze(Bak(state='A', s1=None, s2=None), Bak(state='M'), Bak(state='A'), effector_auto_rates)

    # Pore formation by effectors
    # ------------------------------------
    pore_max_size = 4
    pore_rates = [[2.040816e-04,  # 1.0e-6/v**2
                   1e-3]] * pore_max_size
    assemble_pore_sequential(Bax(bf=None, state='A'), pore_max_size, pore_rates)
    assemble_pore_sequential(Bak(bf=None, state='A'), pore_max_size, pore_rates)

    # ------------------------------------
    # MOMP Inhibition
    # -----------------------------------------------------
    # tBid and free Bax recruit Bcl-xL(C)
    # ------------------------------------------------------------------------
    bclxl_recruitment_rates = [2.040816e-04, 1e-3, 1]
    catalyze(Bid(state='M'), BclxL(state='C'), BclxL(state='M'), bclxl_recruitment_rates)
    catalyze(Bax(state='A', s1=None, s2=None), BclxL(state='C'), BclxL(state='M'), bclxl_recruitment_rates)

    # Bcl2 inhibitors of Bax, Bak, and Bid
    # ------------------------------------
    simple_bind_table([[                                            Bcl2,         BclxL,  Mcl1],
                       [                                              {}, {'state':'M'},    {}],
                       [Bid, {'state':'M'},                         True,          True,  True],
                       [Bax, {'s1':None, 's2':None, 'state':'A'},   True,          True, False],
                       [Bak, {'s1':None, 's2':None, 'state':'A'},  False,          True,  True]],
                      kd['BID_BAX_BAK_inh'], model)

    # Bcl2 sensitizers 
    # ---------------------------------------------------------------
    simple_bind_table([[                      Bcl2,         BclxL,  Mcl1],
                       [                        {}, {'state':'M'},    {}],
                       [Bad,  {'state':'M'},  True,          True, False],
                       [NOXA, {},            False,         False,  True]], 
                      kd['BCLs_sens'], model)

    # CytC, Smac release
    # ----------------------
    #ringp_transport(Subunit, Source, Dest, min_size, max_size, rates):
    ringp_transport(Bax(bf=None), CytoC(state='M'), CytoC(state='C'), 4, 4, kd['BAX_CYTC']) 
    ringp_transport(Bax(bf=None),  Smac(state='M'),  Smac(state='C'), 4, 4, kd['BAX_SMAC']) 
    ringp_transport(Bak(bf=None), CytoC(state='M'), CytoC(state='C'), 4, 4, kd['BAK_CYTC'])
    ringp_transport(Bak(bf=None),  Smac(state='M'),  Smac(state='C'), 4, 4, kd['BAK_SMAC'])

    # --------------------------------------
    # CytC and Smac activation after release
    # --------------------------------------
    Rule('act_cSmac',  Smac(bf=None, state='C') <> Smac(bf=None, state='A'), kd['SMAC_ACT'][0], kd['SMAC_ACT'][1])
    Rule('act_cCytoC', CytoC(bf=None, state='C') <> CytoC(bf=None, state='A'), kd['CYTOC_ACT'][0], kd['CYTOC_ACT'][1])

