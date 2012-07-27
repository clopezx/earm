from pysb.bng import generate_equations

models = ['lopez_direct', 'lopez_indirect', 'lopez_embedded',
          'albeck_11b', 'albeck_11c', 'albeck_11d', 'albeck_11e', 'albeck_11f',
          'chen_2007_biophys_j', 'chen_2007_febs_direct', 'chen_2007_febs_indirect',
          'cui_direct', 'cui_direct1', 'cui_direct2',
          'howells']

for model in models:
    exec('import earm.mito.%s' % model)
    exec('m = earm.mito.%s.model' % model)
    print("----------------------")
    print("Model: %s" % model)
    print("Number of rules: %d" % len(m.rules))
    generate_equations(m)
    print("Number of odes: %d" % len(m.odes))
