from pysb.bng import generate_equations
from pysb import *

models = ['lopez_embedded', 'lopez_direct', 'lopez_indirect',
          'albeck_11b', 'albeck_11c', 'albeck_11d', 'albeck_11e', 'albeck_11f',
          'chen_2007_biophys_j', 'chen_2007_febs_direct', 'chen_2007_febs_indirect',
          'cui_direct', 'cui_direct1', 'cui_direct2',
          'howells']

# Mito only
print('%-24s %11s %10s %12s %11s %10s %12s' %
        ('Model', 'Mito Rules', 'Mito ODEs', 'Mito Params',
         'Full Rules', 'Full ODEs', 'Full Params'))

for model in models:
    exec('import earm.mito.%s' % model)
    exec('import earm.%s' % model)
    exec('m_mito = earm.mito.%s.model' % model)
    exec('m_full = earm.%s.model' % model)

    # Add tBid initial condition to mito only model
    Bid = m_mito.all_components()['Bid']
    try:
        p = Parameter('tBid_0', 1, _export=False)
        m_mito.initial(Bid(state='T', bf=None), p)
        m_mito.add_component(p)
    except Exception:
        pass # Duplicate initial condition for tBid

    if model in ['chen_2007_febs_direct', 'howells']:
        Bad = m_mito.all_components()['Bad']
        p = Parameter('mBad_0', 1, _export=False)
        m_mito.initial(Bad(state='M', bf=None, serine='U'), p)
        m_mito.add_component(p)

    generate_equations(m_mito)
    generate_equations(m_full)
    print("%-24s %11d %10d %12d %11d %10d %12d" %
            (model, len(m_mito.rules), len(m_mito.odes), len(m_mito.parameters),
             len(m_full.rules), len(m_full.odes), len(m_full.parameters)))

# Full models
