"""Figure from Chen 2007, FEBS Letters."""

from pysb import *
from pysb.integrate import odesolve
from pysb.bng import generate_equations
from pylab import linspace, plot, figure, ion, legend, ylim
import re

from earm.mito.chen_2007_febs_indirect import model as indirect
from earm.mito.chen_2007_febs_direct import model as direct

def figure_2a():
    """Reproduce the dose-response in Figure 2a of Chen et al.

    Despite the fact that the PySB versions of the two models exactly reproduce
    the ODEs given in the supplement of the Chen et al. paper, this function
    does not quantitatively reproduce the plots in Figure 2. The steady-state
    output of both models is slightly higher than shown in the figure: for
    example, the steady-state fraction of oligomerized Bax produced by the
    PySB direct model for 0 input stimulus is 0.85, whereas the intercept shown
    in the plot is clearly lower than this. The same is true for the indirect
    model.

    The reasons for this are not clear, as the parameter values and ODE
    structure have been verified to be identical. There appears to be a
    discrepancy in the way that the dose-response curves are calculated,
    possibly in the timepoint used for sampling the steady-state value.
    """

    f_range = linspace(0, 20, 50)
    t = linspace(0, 50000, 1000)
    ion()
    figure()

    for model in [direct, indirect]:
        ss_Bax4_vals = []

        for f in f_range:
            if model is direct:
                model.parameters['Bid_0'].value = 1 + (f * 1)
                model.parameters['Bad_0'].value = 2 + (f * 2)
            else:
                model.parameters['Bid_0'].value = 3 + (f * 3)

            # Calculate the fraction of oligomerized Bax
            x = odesolve(model, t)
            Bax_frac = (4*x['Bax4_'])/model.parameters['Bax_0'].value
            ss_Bax4_val = Bax_frac[-1]
            ss_Bax4_vals.append(ss_Bax4_val)

        plot(f_range, ss_Bax4_vals, label='DM' if model is direct else 'IM')
        ylim([0.75, 1])

    legend(loc='lower right')

