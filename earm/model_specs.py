from pysb.bng import generate_equations

models = ['lopez_direct', 'lopez_indirect', 'lopez_embedded',
          'albeck_11b', 'albeck_11c', 'albeck_11d', 'albeck_11e', 'albeck_11f',
          'chen_2007_biophys_j', 'chen_2007_febs_direct', 'chen_2007_febs_indirect',
          'cui_direct', 'cui_direct1', 'cui_direct2',
          'howells']

print '%-24s %8s %8s' % ('Model', 'Rules', 'ODEs')
for model in models:
    exec('import earm.mito.%s' % model)
    exec('m = earm.mito.%s.model' % model)
    generate_equations(m)
    print "%-24s %8d %8d" % (model, len(m.rules), len(m.odes))
