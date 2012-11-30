"""
This module contains a number of tests that can be used to verify that
the behavior of the MOMP-only sub-models based on [Albeck2008]_
and re-written in PySB match the original publication. To verify that
the PySB version of the models match the original models (which were written
in MATLAB), timecourse data was generated in MATLAB using the original MATLAB
code, and saved in tab-separated data files (`albeck_11{b-f}.tsv`). The PySB
models are then run using the same input and the output is verified to match
the original data (within integration tolerances).

This model verification procedure is written as a series of unit tests, one
for each sub-model, using the built-in Python package unittest.

To run the tests, simply execute this file at the command line, i.e.::

    python test_albeck_models.py

and all tests should pass.
"""

import unittest

from pysb import *
from pysb.bng import generate_equations
from pysb.integrate import odesolve

from earm.mito import albeck_11b
from earm.mito import albeck_11c
from earm.mito import albeck_11d
from earm.mito import albeck_11e
from earm.mito import albeck_11f

from matplotlib.pyplot import figure, ion, plot, legend
import numpy as np
import os.path

# The default integration tolerance in pysb.integrate
rtol = 1e-6

def add_caspase8(model):
    """Add the reaction C8 + Bid <-> C8:Bid -> C8 + tBid.

    All of the MOMP sub-models in this model repository have been written
    to have tBid as their upstream "interface," and Smac/Cytochrome C release
    as their downstream "interface." However, in the original publication
    [Albeck2008]_, the authors incorporated the additional upstream
    step of Bid cleavage by caspase 8 into all of the sub-models that they
    explored.

    Therefore, to match the output of the PySB models in this repository to
    the output produced by the original MATLAB models, it is necessary to add
    the upstream caspase-8:Bid reactions.

    This function takes a model object and adds the necessary elements:
    - Caspase-8 Monomer 
    - Caspase-8 initial condition
    - Caspase-8/Bid cleavage reaction and associated parameters

    In addition, since in the original publication the plotted figures
    considered Smac release kinetics in the absence of Cytochrome C release,
    this function sets the Cytochrome C initial condition to 0 (this prevents
    Cytochrome C from competing with Smac for pore transport, which affects
    the observed Smac release kinetics slightly).
    """

    # If not already added, add upstream caspase reactions to model
    if model.monomers.get('C8') is None:
        Bid = model.monomers.get('Bid')
        # Add caspase 8
        C8 = Monomer('C8', ['state', 'bf'], {'state': ['pro', 'A']},
                     _export=False)
        model.add_component(C8)

        # Add caspase 8 initial condition (with placeholder value 1)
        C8_0 = Parameter('C8_0', 1, _export=False)
        model.add_component(C8_0)
        model.initial(C8(state='A', bf=None), C8_0)

        # Add rules C8 + Bid <-> C8:Bid -> Bid + C8*
        kf = Parameter('bind_C8A_BidU_to_C8ABidU_kf', 1e-7, _export=False)
        kr = Parameter('bind_C8A_BidU_to_C8ABidU_kr', 1e-3, _export=False)
        kc = Parameter('catalyze_C8ABidU_to_C8A_BidT_kc', 1, _export=False)
        model.add_component(kf)
        model.add_component(kr)
        model.add_component(kc)

        rb = Rule('bind_C8A_BidU_to_C8ABidU',
             C8(state='A', bf=None) + Bid(state='U', bf=None) <>
             C8(state='A', bf=1) % Bid(state='U', bf=1),
             kf, kr, _export=False)
        rc = Rule('catalyze_C8ABidU_to_C8A_BidT',
             C8(state='A', bf=1) % Bid(state='U', bf=1) >>
             C8(state='A', bf=None) + Bid(state='T', bf=None),
             kc, _export=False)
        model.add_component(rb)
        model.add_component(rc)

    # Set CytoC to 0 so transport is only of Smac
    model.parameters['CytoC_0'].value = 0

def run_figure_sim(model):
    """Run the C8 dose-response series shown in Fig. 11 of [Albeck2008]_.

    Returns
    -------
    [t, outputs] : list containing two numpy.array objects
        t: The time coordinates of each timepoint, in seconds, from 0 to
            60,000 seconds (15 hours), at 60-second intervals.
        outputs: A 901 x 7 array. Each row corresponds to the fraction of Smac
            released into the cytosol at each timepoint. Each column correspond
            to a distinct caspase 8 dose, from lowest to highest:
            [1, 5, 10, 50, 100, 500, 1000], in molecules per cell.
    """
    # Set timepoints and c8 doses
    tf = 15 * 3600 # 15 hours
    t = np.linspace(0, tf, (tf / 60) + 1)
    c8_doses = [0.01e2, 0.05e2, 0.1e2, 0.5e2, 1e2, 5e2, 10e2]

    outputs = np.empty((len(t), len(c8_doses)))
    for i, c8_dose in enumerate(c8_doses):
        model.parameters['C8_0'].value = c8_dose
        x = odesolve(model, t, rtol=rtol)
        frac_Smac_release = x['cSmac'] / model.parameters['Smac_0'].value
        outputs[:,i] = frac_Smac_release

    return [t, outputs]

def plot_figure(model, data_file):
    """Plot the PySB model output alongside the original MATLAB output.

    This function is not used explicitly by any of the testing code,
    but it is useful for a visual comparison of the output of the PySB
    model to the original MATLAB model output.

    Parameters
    ----------
    model : pysb.model
        The PySB MOMP model.
    data_file : The original MATLAB data file.
    """

    [t, pysb_data] = run_figure_sim(model)
    figure()
    ion()
    pysb_num_doses = pysb_data.shape[1]

    full_data_file_path = os.path.join(os.path.dirname(__file__), data_file)
    mat_data = np.loadtxt(full_data_file_path)
    mat_num_doses = mat_data.shape[1]
    if pysb_num_doses != mat_num_doses:
        raise Exception('Different number of doses between PySB and MATLAB ' +
                        'files.')

    for i in range(pysb_num_doses):
        plot(t, pysb_data[:,i], 'r', label='PySB')
        plot(t, mat_data[:,i], 'b', label='Matlab')

def matches_figure(model, data_file):
    """Test whether the PySB model output matches the original MATLAB output.

    Calls :py:func:`run_figure_sim` to generate the PySB model output, then
    loads the MATLAB data file and compares the outputs using the function
    numpy.allclose.  Returns True if every timepoint from the dose-response
    series matches the original MATLAB output to within one order of magnitude
    of the integration tolerance.
    """

    [t, pysb_data] = run_figure_sim(model)
    full_data_file_path = os.path.join(os.path.dirname(__file__), data_file)
    mat_data = np.loadtxt(full_data_file_path)

    if pysb_data.shape[1] != mat_data.shape[1]:
        raise Exception('Different number of doses between PySB and MATLAB ' +
                        'files.')

    return bool(np.allclose(pysb_data, mat_data, atol=rtol*10))

## TESTS ===============================================================
class TestAlbeck11b(unittest.TestCase):
    """Test the PySB model based on the topology shown in Figure 11b."""
    def setUp(self):
        self.model = albeck_11b.model
        add_caspase8(self.model)

    def test_figure(self):
        """Test that the model reproduces the plots shown in Figure 11b."""
        self.assertTrue(matches_figure(self.model, 'albeck_11b.tsv'))

class TestAlbeck11c(unittest.TestCase):
    """Test the PySB model based on the topology shown in Figure 11c."""
    def setUp(self):
        self.model = albeck_11c.model
        add_caspase8(self.model)

    def test_figure(self):
        """Test that the model reproduces the plots shown in Figure 11c."""
        self.assertTrue(matches_figure(self.model, 'albeck_11c.tsv'))

class TestAlbeck11d(unittest.TestCase):
    """Test the PySB model based on the topology shown in Figure 11d."""
    def setUp(self):
        self.model = albeck_11d.model
        add_caspase8(self.model)

    def test_figure(self):
        """Test that the model reproduces the plots shown in Figure 11d."""
        self.assertTrue(matches_figure(self.model, 'albeck_11d.tsv'))

class TestAlbeck11e(unittest.TestCase):
    """Test the PySB model based on the topology shown in Figure 11e."""
    def setUp(self):
        self.model = albeck_11e.model
        add_caspase8(self.model)

    def test_figure(self):
        """Test that the model reproduces the plots shown in Figure 11e."""
        self.assertTrue(matches_figure(self.model, 'albeck_11e.tsv'))

class TestAlbeck11f(unittest.TestCase):
    """Test the PySB model based on the topology shown in Figure 11f."""
    def setUp(self):
        self.model = albeck_11f.model
        add_caspase8(self.model)

    def test_figure(self):
        """Test that the model reproduces the plots shown in Figure 11f."""
        self.assertTrue(matches_figure(self.model, 'albeck_11f.tsv'))

if __name__ == '__main__':
    unittest.main()
