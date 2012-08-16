# Embedded Model
def embedded():
    """ Direct and indirect modes of action, occurring at the membrane.
    """
    alias_model_components()

    declare_initial_conditions()

    translocate_tBid_Bax_BclxL()

    tBid_activates_Bax_and_Bak()

    tBid_binds_all_anti_apoptotics()

    # Autoactivation: Bax and Bak activate their own kind, but only when
    # free (i.e. not part of a pore complex)
    effector_auto_rates = [1e-7, 1e-3, 1]
    catalyze(Bax(active_monomer), Bax(state='C'), Bax(state='A'),
             effector_auto_rates)
    catalyze(Bak(active_monomer), Bak(state='M'), Bak(state='A'),
             effector_auto_rates)

    # tBid and free Bax recruit Bcl-xL(C)
    bclxl_recruitment_rates = [2.040816e-04, 1e-3, 1]
    catalyze(Bid(state='M'), BclxL(state='C'), BclxL(state='M'),
             bclxl_recruitment_rates)
    catalyze(Bax(active_monomer), BclxL(state='C'), BclxL(state='M'),
             bclxl_recruitment_rates)

    # Anti-apoptotics bind activated effectors
    # Doug Green's "MODE 2" inhibition
    bind_table([[                            Bcl2,  BclxL(state='M'),        Mcl1],
                [Bax(active_monomer),  bcl2_rates,        bcl2_rates,        None],
                [Bak(active_monomer),        None,        bcl2_rates,  bcl2_rates]])

    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition and transport CytoC/Smac
    lopez_pore_formation()

