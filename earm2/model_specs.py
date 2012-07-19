from pysb.bng import generate_equations

models = ['earm2_direct', 'earm2_indirect', 'earm2_embedded',
          'albeck_11b', 'albeck_11c', 'albeck_11d', 'albeck_11e', 'albeck_11f',
          'chen2007BiophysJ', 'chen2007FEBS_direct', 'chen2007FEBS_indirect',
          'cui2008_direct', 'cui2008_direct1', 'cui2008_direct2',
          'howells2011']

for model in models:
    exec('import earm2.mito.%s' % model)
    exec('m = earm2.mito.%s.model' % model)
    print("----------------------")
    print("Model: %s" % model)
    print("Number of rules: %d" % len(m.rules))
    generate_equations(m)
    print("Number of odes: %d" % len(m.odes))
