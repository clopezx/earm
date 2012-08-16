# Embedded Model
def embedded():
    """ Direct and indirect modes of action, occurring at the membrane.
    """
    alias_model_components()
    declare_initial_conditions()

    effector_auto_rates = [1e-7, 1e-3, 1]
    bclxl_recruitment_rates = [2.040816e-04, 1e-3, 1]
    # bid_effector_rates
    # transloc_rates
    # bcl2_rates
    pore_rates = [[2.040816e-04, 1e-3]] * (pore_max_size - 1)

    # tBid, Bax, and BclXL translocate to the membrane
    free_Bax = Bax(bf=None, s1=None, s2=None) # Alias

    equilibrate(Bid(bf=None, state='T'), Bid(bf=None, state='M'), [1e-1, 1e-3])
    equilibrate(free_Bax(state='C'), free_Bax(state='M'), transloc_rates)
    equilibrate(BclxL(bf=None, state='C'), BclxL(bf=None, state='M'), transloc_rates)
    catalyze(Bid(state='M'), Bax(state='M'), Bax(state='A'), bid_effector_rates)
    catalyze(Bid(state='M'), Bak(state='M'), Bak(state='A'), bid_effector_rates)
    bind_table([[                           Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bid(state='M'),      bcl2_rates,        bcl2_rates,       bcl2_rates],
                [Bad(state='M'),      bcl2_rates,        bcl2_rates,             None],
                [Noxa(state='M'),           None,              None,       bcl2_rates],
                [Bax(active_monomer), bcl2_rates,        bcl2_rates,             None],
                [Bak(active_monomer),       None,        bcl2_rates,       bcl2_rates]])
    catalyze(Bax(active_monomer), Bax(state='C'), Bax(state='A'), effector_auto_rates)
    catalyze(Bak(active_monomer), Bak(state='M'), Bak(state='A'), effector_auto_rates)
    catalyze(Bid(state='M'), BclxL(state='C'), BclxL(state='M'),
             bclxl_recruitment_rates)
    catalyze(Bax(active_monomer), BclxL(state='C'), BclxL(state='M'),
             bclxl_recruitment_rates)
    assemble_pore_sequential(Bax(bf=None, state='A'), 4, pore_rates)
    assemble_pore_sequential(Bak(bf=None, state='A'), 4, pore_rates)

    # 1. tBid activates Bax and Bak
    # 2. Binding of tBid, Bad, Noxa, Bax and Bak to Bcl2, Mcl1, and Bcl-XL
    # 3. Autoactivation: Bax and Bak activate their own kind, but only when
    # 4. free (i.e. not part of a pore complex)
    # 5. tBid and free Bax recruit Bcl-xL
    # 6. Bax and Bak form pores by sequential addition


