import unittest
from earm2.mito.chen2007BiophysJ import model as chen2007BiophysJ_model
from earm2.mito.chen2007FEBS_direct import model as chen2007FEBS_direct_model
from earm2.mito.chen2007FEBS_indirect import model \
                                          as chen2007FEBS_indirect_model
from earm2.mito.cui2008_direct import model as cui2008_direct_model
from earm2.mito.cui2008_direct1 import model as cui2008_direct1_model
from earm2.mito.cui2008_direct2 import model as cui2008_direct2_model
from earm2.mito.howells2011 import model as howells2011_model
from pysb.bng import generate_equations
import re

def convert_odes(model, p_name_map, s_name_map):
    name_map = {}
    name_map.update(p_name_map)
    name_map.update(s_name_map)
    generate_equations(model)
    ode_list = []
    for i, ode in enumerate(model.odes):
        new_ode = ode.subs(name_map)
        new_ode = 'd[%s]/dt = %s' % (s_name_map['s%d' % i], str(new_ode))
        ode_list.append(new_ode)

    return ode_list

def chen2007BiophysJ_convert_odes(model):
    """Substitute species/parameter names with ones from Chen et al., Biophys J.

    Parameters
    ----------
    model : pysb.core.Model
        The model derived from Chen et al. (2007) Biophysical Journal.

    Returns
    -------
    A list of strings, with one entry for each ODE in the model. Each ODE
    is represented as a string, e.g. "d[Act]/dt = ..."
    """

    p_name_map = {
        'one_step_BidT_BaxC_to_BidT_BaxA_kf': 'k1',
        'reverse_BaxA_to_BaxC_k': 'k2',
        'bind_BidT_Bcl2_kf': 'k5',
        'bind_BidT_Bcl2_kr': 'k6',
        'bind_BaxA_Bcl2_kf': 'k3',
        'bind_BaxA_Bcl2_kr': 'k4',
        'displace_BaxA_BidTBcl2_to_BaxABcl2_BidT_k': 'k7',
        #'displace_BaxA_BidTBcl2_to_BaxABcl2_BidT_kr': 'k8',
        'spontaneous_pore_BaxA_to_Bax4_kf': 'k9',
        'spontaneous_pore_BaxA_to_Bax4_kr': 'k10' }
    s_name_map = {
        's0': 'Act',
        's1': 'InBax',
        's2': 'Bcl2',
        's3': 'AcBax',
        's4': 'ActBcl2',
        's5': 'AcBaxBcl2',
        's6': 'Bax4'}
    return convert_odes(model, p_name_map, s_name_map)

def chen2007FEBS_convert_odes(model, model_type):
    """Substitute species/parameter names with ones from Chen et al., FEBS.

    Parameters
    ----------
    model : pysb.core.Model
        One of the two models from Chen et al., 2008: "direct" or
        "direct".
    model_type : string
        String specifying which type of model it is: 'direct' or 'indirect'.

    Returns
    -------
    A list of strings, with one entry for each ODE in the model. Each ODE
    is represented as a string, e.g. "d[Act]/dt = ..."
    """

    p_name_map = {
        'one_step_BidT_BaxC_to_BidT_BaxA_kf': 'k_InBax',
        'reverse_BaxA_to_BaxC_k': 'k_Bax',
        'bind_BidT_Bcl2_kf': 'k_BH3_Bcl2',
        'bind_BidT_Bcl2_kr': 'kr_BH3Bcl2',
        'bind_BadM_Bcl2_kf': 'k_BH3_Bcl2',
        'bind_BadM_Bcl2_kr': 'kr_BH3Bcl2',
        'bind_BaxC_Bcl2_kf': 'k_Bax_Bcl2',
        'bind_BaxC_Bcl2_kr': 'kr_BaxBcl2',
        'spontaneous_pore_BaxC_to_Bax4_kf': 'k_o',
        'spontaneous_pore_BaxC_to_Bax4_kr': 'kr_o',
        'spontaneous_pore_BaxA_to_Bax4_kf': 'k_o',
        'spontaneous_pore_BaxA_to_Bax4_kr': 'kr_o' }

    if model_type == 'direct':
        s_name_map = {
            's0': 'Act',
            's1': 'Ena',
            's2': 'InBax',
            's3': 'Bcl2',
            's4': 'Bax',
            's5': 'ActBcl2',
            's6': 'EnaBcl2',
            's7': 'MAC'}
    elif model_type == 'indirect':
        s_name_map = {
            's0': 'BH3',
            's1': 'Bax',
            's2': 'Bcl2',
            's3': 'BH3Bcl2',
            's4': 'BaxBcl2',
            's5': 'MAC'}
    else:
        raise ValueError("model_type must be either 'direct' or 'indirect'.")

    return convert_odes(model, p_name_map, s_name_map)

def cui_convert_odes(model):
    """Substitutes species and parameter names with ones from Cui et al., 2008.

    Parameters
    ----------
    model : pysb.core.Model
        One of the three models from Cui et al., 2008: "direct", "direct 1",
        or "direct 2".

    Returns
    -------
    A list of strings, with one entry for each ODE in the model. Each ODE
    is represented as a string, e.g. "d[Act]/dt = ..."
    """

    # Mapping of parameter names
    p_name_map = {
        'one_step_BidT_BaxC_to_BidT_BaxA_kf': 'k1',
        'reverse_BaxA_to_BaxC_k': 'k8',
        'bind_BidT_Bcl2_kf': 'k4',
        'bind_BidT_Bcl2_kr': 'k5',
        'bind_BadM_Bcl2_kf': 'k9',
        'bind_BadM_Bcl2_kr': 'k10',
        'displace_BaxA_BidTBcl2_to_BaxABcl2_BidT_kf': 'k6',
        'displace_BaxA_BidTBcl2_to_BaxABcl2_BidT_kr': 'k7',
        'displace_BadM_BidTBcl2_to_BadMBcl2_BidT_kf': 'k11',
        'displace_BadM_BidTBcl2_to_BadMBcl2_BidT_kr': 'k12',
        'displace_BadM_BaxABcl2_to_BadMBcl2_BaxA_kf': 'k13',
        'displace_BadM_BaxABcl2_to_BadMBcl2_BaxA_kr': 'k14',
        'dimerize_Bax_kf': 'k16',
        'dimerize_Bax_kr': 'k17',
        'synthesize_BaxC_k': 'p1',
        'degrade_BaxC_k': 'u1',
        'degrade_BaxA_k': 'u2',
        'synthesize_BidT_k': 'p2',
        'degrade_BidT_k': 'u3',
        'synthesize_Bcl2_k': 'p3',
        'degrade_Bcl2_k': 'u4',
        'degrade_BidTBcl2_k': 'u5',
        'degrade_BaxBcl2_k': 'u6',
        'synthesize_BadMU_k': 'p4',
        'degrade_BadMU_k': 'u7',
        'degrade_BadBcl2_k': 'u8',
        'degrade_BaxBax_k': 'u9',
        'bind_BaxA_Bcl2_kf': 'k2',
        'bind_BaxA_Bcl2_kr': 'k3',
        'Bax_autoactivation_dimerization_k': 'k15' }
    # Mapping of species names
    s_name_map = {
        's0': 'Act',
        's1': 'Ena',
        's2': 'InBax',
        's3': 'Bcl2',
        's4': '__source',
        's5': 'AcBax',
        's6': 'ActBcl2',
        's7': 'EnaBcl2',
        's8': '__sink',
        's9': 'MAC',
        's10': 'AcBaxBcl2'}
    return convert_odes(model, p_name_map, s_name_map)

def howells_convert_odes(model):
    """Substitutes species and parameter names with ones from Howells 2011.

    Parameters
    ----------
    model : pysb.core.Model
        The model derived from Howells et al. (2011) Journal of Theoretical
        Biology.

    Returns
    -------
    A list of strings, with one entry for each ODE in the model. Each ODE
    is represented as a string, e.g. "d[Act]/dt = ..."
    """

    # Mapping of parameter names
    p_name_map = {
        'one_step_BidT_BaxC_to_BidT_BaxA_kf': 'k_Bak_cat',
        'reverse_BaxA_to_BaxC_k': 'k_Bak_inac',
        'bind_BidT_Bcl2_kf': 'ka_tBid_Bcl2',
        'bind_BidT_Bcl2_kr': 'kd_tBid_Bcl2',
        'bind_BaxA_Bcl2_kf': 'ka_Bak_Bcl2',
        'bind_BaxA_Bcl2_kr': 'kd_Bak_Bcl2',
        'displace_BaxA_BidTBcl2_to_BaxABcl2_BidT_k': 'k_tBid_rel2',
        'spontaneous_pore_BaxA_to_Bax4_kf': 'ka_Bak_poly',
        'spontaneous_pore_BaxA_to_Bax4_kr': 'kd_Bak_poly',
        'equilibrate_BadCU_to_BadMU_kf': 't_Bad_in',
        'equilibrate_BadCU_to_BadMU_kr': 't_Bad_out',
        'bind_BadM_Bcl2_kf': 'ka_Bad_Bcl2',
        'bind_BadM_Bcl2_kr': 'kd_Bad_Bcl2',
        'displace_BadM_BidTBcl2_to_BadMBcl2_BidT_k': 'k_tBid_rel1',
        'phosphorylate_BadC_k': 'k_Bad_phos1',
        'phosphorylate_BadM_k': 'k_Bad_phos1',
        'phosphorylate_BadMUBcl2_to_BadCP_k': 'k_Bad_phos2',
        'sequester_BadCP_to_BadC1433_k': 'k_Bad_seq',
        'release_BadC1433_to_BadCU_k': 'k_Bad_rel'}
    # Mapping of species names
    s_name_map = {
        's0': 'tBid',
        's1': 'Bak_inac',
        's2': 'Bcl2',
        's3': 'Bad_m',
        's4': 'Bak',
        's5': 'tBidBcl2',
        's6': 'Bad',
        's7': 'BadBcl2',
        's8': 'pBad',
        's9': 'BakBcl2',
        's10': 'Bak_poly',
        's11': 'pBad1433'}
    return convert_odes(model, p_name_map, s_name_map)

## TESTS ===============================================================

class TestChen2007BiophysJ(unittest.TestCase):
    def setUp(self):
        self.model = chen2007BiophysJ_model

    def test_odes(self):
        """These ODEs match those of the paper with the caveat that in the
        equation for d[Bcl2]/dt, in the paper the authors group the terms

            + AcBaxBcl2*k4 + ActBcl2*k6'

        into the single term

            + k_Bcl2 * Bcl2_nonfree

        with the comment that "Bcl2_nonfree indicates the total concentration
        of Bcl2 associated with both activated Bax and Activator
        ([Bcl2_nonfree] = [AcBaxBcl2] + [ActBcl2]). We use k_bcl2 to represent
        the rate of non-free Bcl2 shifting to free Bcl2, assuming that free
        Bcl2 originates from both Bcl2 non-free forms at the same rate."

        In addition, in the legend for Table 1 (which lists parameters) they
        state that: "We set k_bcl2 (the rate of non-free Bcl2 shifting to
        free Bcl2)... equal to k6 assuming that free Bcl2 originate [sic] from
        AcBaxBcl2 at the same rate with ActBcl2."

        So, obviously this substitution of Bcl2_nonfree for AcBaxBcl2 and
        ActBcl2 works if k4 and k6 are equal, which they claim as an
        assumption; however, in their table of parameter values, they list k4
        as having a value of 0.001 s^-1, and k6 as having a value of 0.04 s^-1.
        """

        ode_list = chen2007BiophysJ_convert_odes(self.model)
        self.assertEqual(ode_list,
            ['d[Act]/dt = AcBax*ActBcl2*k7 - Act*Bcl2*k5 + ActBcl2*k6',
             'd[InBax]/dt = AcBax*k2 - Act*InBax*k1',
             'd[Bcl2]/dt = -AcBax*Bcl2*k3 + AcBaxBcl2*k4 - Act*Bcl2*k5 + ActBcl2*k6',
             'd[AcBax]/dt = -1.0*AcBax**4*k9 - AcBax*ActBcl2*k7 - AcBax*Bcl2*k3 - AcBax*k2 + AcBaxBcl2*k4 + Act*InBax*k1 + 4*Bax4*k10',
             'd[ActBcl2]/dt = -AcBax*ActBcl2*k7 + Act*Bcl2*k5 - ActBcl2*k6',
             'd[AcBaxBcl2]/dt = AcBax*ActBcl2*k7 + AcBax*Bcl2*k3 - AcBaxBcl2*k4',
             'd[Bax4]/dt = 0.25*AcBax**4*k9 - Bax4*k10'])

class TestChen2007FEBS_Indirect(unittest.TestCase):
    def setUp(self):
        self.model = chen2007FEBS_indirect_model

    def test_odes(self):
        ode_list = chen2007FEBS_convert_odes(self.model, 'indirect')
        self.assertEqual(ode_list,
            ['d[BH3]/dt = -BH3*Bcl2*k_BH3_Bcl2 + BH3Bcl2*kr_BH3Bcl2',
             'd[Bax]/dt = -1.0*Bax**4*k_o - Bax*Bcl2*k_Bax_Bcl2 + BaxBcl2*kr_BaxBcl2 + 4*MAC*kr_o',
             'd[Bcl2]/dt = -BH3*Bcl2*k_BH3_Bcl2 + BH3Bcl2*kr_BH3Bcl2 - Bax*Bcl2*k_Bax_Bcl2 + BaxBcl2*kr_BaxBcl2',
             'd[BH3Bcl2]/dt = BH3*Bcl2*k_BH3_Bcl2 - BH3Bcl2*kr_BH3Bcl2',
             'd[BaxBcl2]/dt = Bax*Bcl2*k_Bax_Bcl2 - BaxBcl2*kr_BaxBcl2',
             'd[MAC]/dt = 0.25*Bax**4*k_o - MAC*kr_o'])

class TestChen2007FEBS_Direct(unittest.TestCase):
    def setUp(self):
        self.model = chen2007FEBS_direct_model

    def test_odes(self):
        ode_list = chen2007FEBS_convert_odes(self.model, 'direct')
        self.assertEqual(ode_list,
            ['d[Act]/dt = -Act*Bcl2*k_BH3_Bcl2 + ActBcl2*kr_BH3Bcl2',
             'd[Ena]/dt = -Bcl2*Ena*k_BH3_Bcl2 + EnaBcl2*kr_BH3Bcl2',
             'd[InBax]/dt = -Act*InBax*k_InBax + Bax*k_Bax',
             'd[Bcl2]/dt = -Act*Bcl2*k_BH3_Bcl2 + ActBcl2*kr_BH3Bcl2 - Bcl2*Ena*k_BH3_Bcl2 + EnaBcl2*kr_BH3Bcl2',
             'd[Bax]/dt = Act*InBax*k_InBax - 1.0*Bax**4*k_o - Bax*k_Bax + 4*MAC*kr_o',
             'd[ActBcl2]/dt = Act*Bcl2*k_BH3_Bcl2 - ActBcl2*kr_BH3Bcl2',
             'd[EnaBcl2]/dt = Bcl2*Ena*k_BH3_Bcl2 - EnaBcl2*kr_BH3Bcl2',
             'd[MAC]/dt = 0.25*Bax**4*k_o - MAC*kr_o'])

class TestCui2008_Direct(unittest.TestCase):
    def setUp(self):
        self.model = cui2008_direct_model

    def test_odes(self):
        ode_list = cui_convert_odes(self.model)
        self.assertEqual(ode_list,
            ['d[Act]/dt = -Act*Bcl2*k4 - Act*EnaBcl2*k12 - Act*u3 + ActBcl2*Ena*k11 + ActBcl2*k5 + __source*p2',
             'd[Ena]/dt = Act*EnaBcl2*k12 - ActBcl2*Ena*k11 - Bcl2*Ena*k9 - Ena*u7 + EnaBcl2*k10 + __source*p4',
             'd[InBax]/dt = AcBax*k8 - Act*InBax*k1 - InBax*u1 + __source*p1',
             'd[Bcl2]/dt = -Act*Bcl2*k4 + ActBcl2*k5 - Bcl2*Ena*k9 - Bcl2*u4 + EnaBcl2*k10 + __source*p3',
             'd[__source]/dt = 0',
             'd[AcBax]/dt = -1.0*AcBax**2*k16 - AcBax*k8 - AcBax*u2 + Act*InBax*k1 + 2*MAC*k17',
             'd[ActBcl2]/dt = Act*Bcl2*k4 + Act*EnaBcl2*k12 - ActBcl2*Ena*k11 - ActBcl2*k5 - ActBcl2*u5',
             'd[EnaBcl2]/dt = -Act*EnaBcl2*k12 + ActBcl2*Ena*k11 + Bcl2*Ena*k9 - EnaBcl2*k10 - EnaBcl2*u8',
             'd[__sink]/dt = AcBax*u2 + Act*u3 + ActBcl2*u5 + Bcl2*u4 + Ena*u7 + EnaBcl2*u8 + InBax*u1 + MAC*u9',
             'd[MAC]/dt = 0.5*AcBax**2*k16 - MAC*k17 - MAC*u9'])

class TestCui2008_Direct1(unittest.TestCase):
    def setUp(self):
        self.model = cui2008_direct1_model

    def test_odes(self):
        ode_list = cui_convert_odes(self.model)
        self.assertEqual(ode_list,
            ['d[Act]/dt = AcBax*ActBcl2*k6 - AcBaxBcl2*Act*k7 - Act*Bcl2*k4 - Act*EnaBcl2*k12 - Act*u3 + ActBcl2*Ena*k11 + ActBcl2*k5 + __source*p2',
             'd[Ena]/dt = AcBax*EnaBcl2*k14 - AcBaxBcl2*Ena*k13 + Act*EnaBcl2*k12 - ActBcl2*Ena*k11 - Bcl2*Ena*k9 - Ena*u7 + EnaBcl2*k10 + __source*p4',
             'd[InBax]/dt = AcBax*k8 - Act*InBax*k1 - InBax*u1 + __source*p1',
             'd[Bcl2]/dt = -AcBax*Bcl2*k2 + AcBaxBcl2*k3 - Act*Bcl2*k4 + ActBcl2*k5 - Bcl2*Ena*k9 - Bcl2*u4 + EnaBcl2*k10 + __source*p3',
             'd[__source]/dt = 0',
             'd[AcBax]/dt = -1.0*AcBax**2*k16 - AcBax*ActBcl2*k6 - AcBax*Bcl2*k2 - AcBax*EnaBcl2*k14 - AcBax*k8 - AcBax*u2 + AcBaxBcl2*Act*k7 + AcBaxBcl2*Ena*k13 + AcBaxBcl2*k3 + Act*InBax*k1 + 2*MAC*k17',
             'd[ActBcl2]/dt = -AcBax*ActBcl2*k6 + AcBaxBcl2*Act*k7 + Act*Bcl2*k4 + Act*EnaBcl2*k12 - ActBcl2*Ena*k11 - ActBcl2*k5 - ActBcl2*u5',
             'd[EnaBcl2]/dt = -AcBax*EnaBcl2*k14 + AcBaxBcl2*Ena*k13 - Act*EnaBcl2*k12 + ActBcl2*Ena*k11 + Bcl2*Ena*k9 - EnaBcl2*k10 - EnaBcl2*u8',
             'd[__sink]/dt = AcBax*u2 + AcBaxBcl2*u6 + Act*u3 + ActBcl2*u5 + Bcl2*u4 + Ena*u7 + EnaBcl2*u8 + InBax*u1 + MAC*u9',
             'd[MAC]/dt = 0.5*AcBax**2*k16 - MAC*k17 - MAC*u9',
             'd[AcBaxBcl2]/dt = AcBax*ActBcl2*k6 + AcBax*Bcl2*k2 + AcBax*EnaBcl2*k14 - AcBaxBcl2*Act*k7 - AcBaxBcl2*Ena*k13 - AcBaxBcl2*k3 - AcBaxBcl2*u6'])

class TestCui2008_Direct2(unittest.TestCase):
    def setUp(self):
        self.model = cui2008_direct2_model

    def test_odes(self):
        ode_list = cui_convert_odes(self.model)
        self.assertEqual(ode_list,
            ['d[Act]/dt = AcBax*ActBcl2*k6 - AcBaxBcl2*Act*k7 - Act*Bcl2*k4 - Act*EnaBcl2*k12 - Act*u3 + ActBcl2*Ena*k11 + ActBcl2*k5 + __source*p2',
             'd[Ena]/dt = AcBax*EnaBcl2*k14 - AcBaxBcl2*Ena*k13 + Act*EnaBcl2*k12 - ActBcl2*Ena*k11 - Bcl2*Ena*k9 - Ena*u7 + EnaBcl2*k10 + __source*p4',
             'd[InBax]/dt = -AcBax*InBax*k15 + AcBax*k8 - Act*InBax*k1 - InBax*u1 + __source*p1',
             'd[Bcl2]/dt = -AcBax*Bcl2*k2 + AcBaxBcl2*k3 - Act*Bcl2*k4 + ActBcl2*k5 - Bcl2*Ena*k9 - Bcl2*u4 + EnaBcl2*k10 + __source*p3',
             'd[__source]/dt = 0',
             'd[AcBax]/dt = -1.0*AcBax**2*k16 - AcBax*ActBcl2*k6 - AcBax*Bcl2*k2 - AcBax*EnaBcl2*k14 - AcBax*InBax*k15 - AcBax*k8 - AcBax*u2 + AcBaxBcl2*Act*k7 + AcBaxBcl2*Ena*k13 + AcBaxBcl2*k3 + Act*InBax*k1 + 2*MAC*k17',
             'd[ActBcl2]/dt = -AcBax*ActBcl2*k6 + AcBaxBcl2*Act*k7 + Act*Bcl2*k4 + Act*EnaBcl2*k12 - ActBcl2*Ena*k11 - ActBcl2*k5 - ActBcl2*u5',
             'd[EnaBcl2]/dt = -AcBax*EnaBcl2*k14 + AcBaxBcl2*Ena*k13 - Act*EnaBcl2*k12 + ActBcl2*Ena*k11 + Bcl2*Ena*k9 - EnaBcl2*k10 - EnaBcl2*u8',
             'd[__sink]/dt = AcBax*u2 + AcBaxBcl2*u6 + Act*u3 + ActBcl2*u5 + Bcl2*u4 + Ena*u7 + EnaBcl2*u8 + InBax*u1 + MAC*u9',
             'd[MAC]/dt = 0.5*AcBax**2*k16 + AcBax*InBax*k15 - MAC*k17 - MAC*u9',
             'd[AcBaxBcl2]/dt = AcBax*ActBcl2*k6 + AcBax*Bcl2*k2 + AcBax*EnaBcl2*k14 - AcBaxBcl2*Act*k7 - AcBaxBcl2*Ena*k13 - AcBaxBcl2*k3 - AcBaxBcl2*u6'])

# TODO: Verify correctness of output
class TestHowells2011(unittest.TestCase):
    def setUp(self):
        self.model = howells2011_model

    def test_odes(self):
        ode_list = howells_convert_odes(self.model)
        self.assertEqual(ode_list,
            ['d[tBid]/dt = Bad_m*k_tBid_rel1*tBidBcl2 + Bak*k_tBid_rel2*tBidBcl2 - Bcl2*ka_tBid_Bcl2*tBid + kd_tBid_Bcl2*tBidBcl2',
             'd[Bak_inac]/dt = Bak*k_Bak_inac - Bak_inac*k_Bak_cat*tBid',
             'd[Bcl2]/dt = BadBcl2*k_Bad_phos2 + BadBcl2*kd_Bad_Bcl2 - Bad_m*Bcl2*ka_Bad_Bcl2 - Bak*Bcl2*ka_Bak_Bcl2 + BakBcl2*kd_Bak_Bcl2 - Bcl2*ka_tBid_Bcl2*tBid + kd_tBid_Bcl2*tBidBcl2',
             'd[Bad_m]/dt = Bad*t_Bad_in + BadBcl2*kd_Bad_Bcl2 - Bad_m*Bcl2*ka_Bad_Bcl2 - Bad_m*k_Bad_phos1 - Bad_m*k_tBid_rel1*tBidBcl2 - Bad_m*t_Bad_out',
             'd[Bak]/dt = -1.0*Bak**4*ka_Bak_poly - Bak*Bcl2*ka_Bak_Bcl2 - Bak*k_Bak_inac - Bak*k_tBid_rel2*tBidBcl2 + BakBcl2*kd_Bak_Bcl2 + Bak_inac*k_Bak_cat*tBid + 4*Bak_poly*kd_Bak_poly',
             'd[tBidBcl2]/dt = -Bad_m*k_tBid_rel1*tBidBcl2 - Bak*k_tBid_rel2*tBidBcl2 + Bcl2*ka_tBid_Bcl2*tBid - kd_tBid_Bcl2*tBidBcl2',
             'd[Bad]/dt = -Bad*k_Bad_phos1 - Bad*t_Bad_in + Bad_m*t_Bad_out + k_Bad_rel*pBad1433',
             'd[BadBcl2]/dt = -BadBcl2*k_Bad_phos2 - BadBcl2*kd_Bad_Bcl2 + Bad_m*Bcl2*ka_Bad_Bcl2 + Bad_m*k_tBid_rel1*tBidBcl2',
             'd[pBad]/dt = Bad*k_Bad_phos1 + BadBcl2*k_Bad_phos2 + Bad_m*k_Bad_phos1 - k_Bad_seq*pBad',
             'd[BakBcl2]/dt = Bak*Bcl2*ka_Bak_Bcl2 + Bak*k_tBid_rel2*tBidBcl2 - BakBcl2*kd_Bak_Bcl2',
             'd[Bak_poly]/dt = 0.25*Bak**4*ka_Bak_poly - Bak_poly*kd_Bak_poly',
             'd[pBad1433]/dt = -k_Bad_rel*pBad1433 + k_Bad_seq*pBad'])


if __name__ == '__main__':
    unittest.main()
