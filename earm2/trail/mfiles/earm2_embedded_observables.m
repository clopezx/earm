function [ode_observables, kd_values, kd_index, ic_index, dividing_factor] = earm2_embedded_observables()

ode_observables{1, 1} = [11 28];
ode_observables{1, 2} = [1 1];
ode_observables{1, 3} = 1;
ode_observables{2, 1} = [31];
ode_observables{2, 2} = [1];
ode_observables{2, 3} = 1;
ode_observables{3, 1} = [68 71];
ode_observables{3, 2} = [1 1];
ode_observables{3, 3} = 1;
ode_observables{4, 1} = [20 63 65];
ode_observables{4, 2} = [1 1 1];
ode_observables{4, 3} = 1;
ode_observables{5, 1} = [38];
ode_observables{5, 2} = [1];
ode_observables{5, 3} = 1;
end
