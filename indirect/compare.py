from earm2.trail.earm2_indirect import model as new_model
from earm_indirect_mem import model as old_model

from pysb.integrate import odesolve
import numpy
import matplotlib.pyplot as plt
from pylab import ion

t = numpy.linspace(0, 8000, 100)
x_new = odesolve(new_model, t)
x_old = odesolve(old_model, t)

print('old indirect rules: %d' % len(old_model.rules))
print('new indirect rules: %d' % len(new_model.rules))
print
print('old indirect species: %d' % len(old_model.species))
print('new indirect species: %d' % len(new_model.species))
print
print('old indirect parameters: %d' % len(old_model.parameters))
print('new indirect parameters: %d' % len(new_model.parameters))
print
print('old indirect reactions: %d' % len(old_model.reactions))
print('new indirect reactions: %d' % len(new_model.reactions))
print
print('old indirect initial conditions: %d' % len(old_model.initial_conditions))
print('new indirect initial conditions: %d' % len(new_model.initial_conditions))

p_old = sorted([(p.value, p.name) for p in old_model.parameters])
p_new = sorted([(p.value, p.name) for p in new_model.parameters])

rates_old = [a['rate'] for a in old_model.reactions]
rates_new = [a['rate'] for a in new_model.reactions]

# Diff will not work b/c of extra parameters
#diff = numpy.array([x[0] for x in p_new]) - numpy.array([x[0] for x in p_old])

#f_old = open('old_params.txt', 'w')
#f_old.write('\n'.join([str(p) for p in p_old]))
#f_old.close()

#f_new = open('new_params.txt', 'w')
#f_new.write('\n'.join([str(p) for p in p_new]))
#f_new.close()

#print 'old model species: %d' % len(old_model.species)
#print 'new model species: %d' % len(new_model.species)

ion()
plt.figure()
plt.plot(t, x_old['cPARP'], 'r')
plt.plot(t, x_new['cPARP'], 'g')



