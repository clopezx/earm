#! /usr/bin/python

from pysb.integrate import odesolve
from pylab import *

from BioPhysJ2007 import model

def run_model(plot_title):
    Act_0 = model.parameter('Act_0').value
    Bax_0 = model.parameter('Bax_0').value
    Bcl2_0 = model.parameter('Bcl2_0').value

    #hrs = 24
    #hrs = 24
    #t = linspace(0, hrs*3600, hrs*60+1)
    t = linspace(0, 3000, 201)

    x = odesolve(model, t)

    #plot(t/3600, x['Act'] / Act_0, label='Act/Act_0', color='r')
    #plot(t/3600, x['AcBax_free'], label='AcBax_free/Bax_0', color='r', linewidth=2)
    plot(t, x['AcBax_tot'] / Bax_0, label='AcBax_tot/Bax_0', color='r', linewidth=2)
    #plot(t, x['InBax_tot'], label='InBax_tot/Bax_0', color='g', linewidth=2)
    #plot(t/3600, x['AcBax_tot'], label='AcBax_tot', color='g', linewidth=2)
    #plot(t/3600, x['AcBax_Bcl2'], label='AcBax_Bcl2/Bax_0', color='b', linewidth=2)
    #plot(t, x['Bcl2_free'], label='Bcl2/Bax_0', color='b', linewidth=2)
    #plot(t/3600, (x['Bax4']*4 / Bax_0), label='Bax4/Bax_0', color='m', linewidth=2)

    xlabel('Time (Hours)')
    ylabel('Species Concentration (Normalized)')
    title(plot_title)
    #legend(loc='NorthEastOutside')
    #ylim(-0.005, 1.005)

run_model('Test')


