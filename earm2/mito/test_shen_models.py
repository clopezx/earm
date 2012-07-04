import unittest
from earm2.mito.cui2008_direct import model as cui2008_direct_model
from earm2.mito.cui2008_direct1 import model as cui2008_direct1_model
from earm2.mito.cui2008_direct2 import model as cui2008_direct2_model
from earm2.mito.chen2007FEBS_direct import model as chen2007FEBS_direct_model
from earm2.mito.chen2007FEBS_indirect import model \
                                          as chen2007FEBS_indirect_model
from pysb.bng import generate_equations
import re

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

    generate_equations(model)
    ode_list = []
    for i, ode in enumerate(model.odes):
        new_ode = 'd[s%d]/dt = %s' % (i, str(ode))
        for old_s in s_name_map:
            new_ode = re.sub(old_s, s_name_map[old_s], new_ode)
        for old_p in p_name_map:
            new_ode = re.sub(old_p, p_name_map[old_p], new_ode)
        ode_list.append(new_ode)

    return ode_list

def cui_convert_odes(model):
    """Substitutes species and parameter names with ones from Cui 2008.

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
        'synthesize_BadM_k': 'p4',
        'degrade_BadM_k': 'u7',
        'degrade_BadBcl2_k': 'u8',
        'degrade_BaxBax_k': 'u9',
        'bind_BaxA_Bcl2_kf': 'k2',
        'bind_BaxA_Bcl2_kr': 'k3',
        'one_step_BaxA_BaxC_to_BaxA_BaxA_kf': 'k15' }
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
        's10': 'MAC',
        's9': 'AcBaxBcl2'}
    generate_equations(model)
    ode_list = []
    for i, ode in enumerate(model.odes):
        new_ode = 'd[s%d]/dt = %s' % (i, str(ode))
        for old_s in s_name_map:
            new_ode = re.sub(old_s, s_name_map[old_s], new_ode)
        for old_p in p_name_map:
            new_ode = re.sub(old_p, p_name_map[old_p], new_ode)
        ode_list.append(new_ode)

    return ode_list

# TODO: Verify correctness of output
class TestChen2007FEBS_Indirect(unittest.TestCase):
    def setUp(self):
        self.model = chen2007FEBS_indirect_model

    def test_odes(self):
        ode_list = chen2007FEBS_convert_odes(self.model, 'indirect')
        self.assertEqual(ode_list,
            ['d[BH3]/dt = -k_BH3_Bcl2*BH3*Bcl2 + kr_BH3Bcl2*BH3Bcl2',
             'd[Bax]/dt = -k_Bax_Bcl2*Bax*Bcl2 + kr_BaxBcl2*BaxBcl2 - 1.0*Bax**4*spontaneous_pore_BaxC_to_Bax4_kf + 4*MAC*spontaneous_pore_BaxC_to_Bax4_kr',
             'd[Bcl2]/dt = -k_Bax_Bcl2*Bax*Bcl2 + kr_BaxBcl2*BaxBcl2 - k_BH3_Bcl2*BH3*Bcl2 + kr_BH3Bcl2*BH3Bcl2',
             'd[BH3Bcl2]/dt = k_BH3_Bcl2*BH3*Bcl2 - kr_BH3Bcl2*BH3Bcl2',
             'd[BaxBcl2]/dt = k_Bax_Bcl2*Bax*Bcl2 - kr_BaxBcl2*BaxBcl2',
             'd[MAC]/dt = 0.25*Bax**4*spontaneous_pore_BaxC_to_Bax4_kf - MAC*spontaneous_pore_BaxC_to_Bax4_kr'])

# TODO: Verify correctness of output
class TestChen2007FEBS_Direct(unittest.TestCase):
    def setUp(self):
        self.model = chen2007FEBS_direct_model

    def test_odes(self):
        ode_list = chen2007FEBS_convert_odes(self.model, 'direct')
        self.assertEqual(ode_list,
            ['d[Act]/dt = -k_BH3_Bcl2*Act*Bcl2 + kr_BH3Bcl2*ActBcl2',
             'd[Ena]/dt = -k_BH3_Bcl2*Ena*Bcl2 + kr_BH3Bcl2*EnaBcl2',
             'd[InBax]/dt = -k_InBax*Act*InBax + k_Bax*Bax',
             'd[Bcl2]/dt = -k_BH3_Bcl2*Ena*Bcl2 + kr_BH3Bcl2*EnaBcl2 - k_BH3_Bcl2*Act*Bcl2 + kr_BH3Bcl2*ActBcl2',
             'd[Bax]/dt = k_InBax*Act*InBax - k_Bax*Bax - 1.0*Bax**4*k_o + 4*MAC*kr_o',
             'd[ActBcl2]/dt = k_BH3_Bcl2*Act*Bcl2 - kr_BH3Bcl2*ActBcl2',
             'd[EnaBcl2]/dt = k_BH3_Bcl2*Ena*Bcl2 - kr_BH3Bcl2*EnaBcl2',
             'd[MAC]/dt = 0.25*Bax**4*k_o - MAC*kr_o'])

# TODO: Verify correctness of output
class TestCui2008_Direct(unittest.TestCase):
    def setUp(self):
        self.model = cui2008_direct_model

    def test_odes(self):
        ode_list = cui_convert_odes(self.model)
        self.assertEqual(ode_list,
            ['d[Act]/dt = -k4*Act*Bcl2 + k5*ActBcl2 - u3*Act + k11*Ena*ActBcl2 - k12*Act*EnaBcl2 + __source*p2',
             'd[Ena]/dt = -k9*Ena*Bcl2 + k10*EnaBcl2 - u7*Ena - k11*Ena*ActBcl2 + k12*Act*EnaBcl2 + __source*p4',
             'd[InBax]/dt = -u1*InBax - k1*Act*InBax + k8*AcBax + __source*p1',
             'd[Bcl2]/dt = -k9*Ena*Bcl2 + k10*EnaBcl2 - k4*Act*Bcl2 + k5*ActBcl2 - u4*Bcl2 + __source*p3',
             'd[__source]/dt = 0',
             'd[AcBax]/dt = -u2*AcBax - 1.0*k16*AcBax**2 + 2*k17*AcBaxBcl2 + k1*Act*InBax - k8*AcBax',
             'd[ActBcl2]/dt = k4*Act*Bcl2 - k5*ActBcl2 - u5*ActBcl2 - k11*Ena*ActBcl2 + k12*Act*EnaBcl2',
             'd[EnaBcl2]/dt = k9*Ena*Bcl2 - k10*EnaBcl2 - u8*EnaBcl2 + k11*Ena*ActBcl2 - k12*Act*EnaBcl2',
             'd[__sink]/dt = u8*EnaBcl2 + u7*Ena + u2*AcBax + u9*AcBaxBcl2 + u1*InBax + u4*Bcl2 + u5*ActBcl2 + u3*Act',
             'd[AcBaxBcl2]/dt = -u9*AcBaxBcl2 + 0.5*k16*AcBax**2 - k17*AcBaxBcl2'])

# TODO: Verify correctness of output
class TestCui2008_Direct1(unittest.TestCase):
    def setUp(self):
        self.model = cui2008_direct1_model

    def test_odes(self):
        ode_list = cui_convert_odes(self.model)
        self.assertEqual(ode_list,
            ['d[Act]/dt = -k4*Act*Bcl2 + k5*ActBcl2 - u3*Act + k11*Ena*ActBcl2 - k12*Act*EnaBcl2 + k6*AcBax*ActBcl2 - k7*Act*MAC + __source*p2',
             'd[Ena]/dt = -k9*Ena*Bcl2 + k10*EnaBcl2 - u7*Ena - k13*Ena*MAC + k14*AcBax*EnaBcl2 - k11*Ena*ActBcl2 + k12*Act*EnaBcl2 + __source*p4',
             'd[InBax]/dt = -u1*InBax - k1*Act*InBax + k8*AcBax + __source*p1',
             'd[Bcl2]/dt = -k9*Ena*Bcl2 + k10*EnaBcl2 - k2*Bcl2*AcBax + k3*MAC - k4*Act*Bcl2 + k5*ActBcl2 - u4*Bcl2 + __source*p3',
             'd[__source]/dt = 0',
             'd[AcBax]/dt = -k2*Bcl2*AcBax + k3*MAC - u2*AcBax - 1.0*k16*AcBax**2 + 2*k17*AcBaxBcl2 + k13*Ena*MAC - k14*AcBax*EnaBcl2 - k6*AcBax*ActBcl2 + k7*Act*MAC + k1*Act*InBax - k8*AcBax',
             'd[ActBcl2]/dt = k4*Act*Bcl2 - k5*ActBcl2 - u5*ActBcl2 - k11*Ena*ActBcl2 + k12*Act*EnaBcl2 - k6*AcBax*ActBcl2 + k7*Act*MAC',
             'd[EnaBcl2]/dt = k9*Ena*Bcl2 - k10*EnaBcl2 - u8*EnaBcl2 + k13*Ena*MAC - k14*AcBax*EnaBcl2 + k11*Ena*ActBcl2 - k12*Act*EnaBcl2',
             'd[__sink]/dt = u8*EnaBcl2 + u7*Ena + u2*AcBax + u9*AcBaxBcl2 + u6*MAC + u1*InBax + u4*Bcl2 + u5*ActBcl2 + u3*Act',
             'd[AcBaxBcl2]/dt = -u9*AcBaxBcl2 + 0.5*k16*AcBax**2 - k17*AcBaxBcl2',
             'd[MAC]/dt = k2*Bcl2*AcBax - k3*MAC - u6*MAC - k13*Ena*MAC + k14*AcBax*EnaBcl2 + k6*AcBax*ActBcl2 - k7*Act*MAC'])

# TODO: Verify correctness of output
class TestCui2008_Direct2(unittest.TestCase):
    def setUp(self):
        self.model = cui2008_direct2_model

    def test_odes(self):
        ode_list = cui_convert_odes(self.model)
        self.assertEqual(ode_list,
            ['d[Act]/dt = -k4*Act*Bcl2 + k5*ActBcl2 - u3*Act + k11*Ena*ActBcl2 - k12*Act*EnaBcl2 + k6*AcBax*ActBcl2 - k7*Act*MAC + __source*p2',
             'd[Ena]/dt = -k9*Ena*Bcl2 + k10*EnaBcl2 - u7*Ena - k13*Ena*MAC + k14*AcBax*EnaBcl2 - k11*Ena*ActBcl2 + k12*Act*EnaBcl2 + __source*p4',
             'd[InBax]/dt = -u1*InBax - k15*InBax*AcBax - k1*Act*InBax + k8*AcBax + __source*p1',
             'd[Bcl2]/dt = -k9*Ena*Bcl2 + k10*EnaBcl2 - k2*Bcl2*AcBax + k3*MAC - k4*Act*Bcl2 + k5*ActBcl2 - u4*Bcl2 + __source*p3',
             'd[__source]/dt = 0',
             'd[AcBax]/dt = -k2*Bcl2*AcBax + k3*MAC - u2*AcBax - 1.0*k16*AcBax**2 + 2*k17*AcBaxBcl2 + k13*Ena*MAC - k14*AcBax*EnaBcl2 - k6*AcBax*ActBcl2 + k7*Act*MAC + k15*InBax*AcBax + k1*Act*InBax - k8*AcBax',
             'd[ActBcl2]/dt = k4*Act*Bcl2 - k5*ActBcl2 - u5*ActBcl2 - k11*Ena*ActBcl2 + k12*Act*EnaBcl2 - k6*AcBax*ActBcl2 + k7*Act*MAC',
             'd[EnaBcl2]/dt = k9*Ena*Bcl2 - k10*EnaBcl2 - u8*EnaBcl2 + k13*Ena*MAC - k14*AcBax*EnaBcl2 + k11*Ena*ActBcl2 - k12*Act*EnaBcl2',
             'd[__sink]/dt = u8*EnaBcl2 + u7*Ena + u2*AcBax + u9*AcBaxBcl2 + u6*MAC + u1*InBax + u4*Bcl2 + u5*ActBcl2 + u3*Act',
             'd[AcBaxBcl2]/dt = -u9*AcBaxBcl2 + 0.5*k16*AcBax**2 - k17*AcBaxBcl2',
             'd[MAC]/dt = k2*Bcl2*AcBax - k3*MAC - u6*MAC - k13*Ena*MAC + k14*AcBax*EnaBcl2 + k6*AcBax*ActBcl2 - k7*Act*MAC'])

if __name__ == '__main__':
    unittest.main()
