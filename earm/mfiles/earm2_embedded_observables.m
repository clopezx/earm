function [ode_observables, kd_values, kd_index, ic_index, dividing_factor] = earm2_embedded_observables()

ode_observables{1, 1} = [12 27]; % Observable('Bid_', Bid(state=U))
ode_observables{1, 2} = [1 1];
ode_observables{1, 3} = 1;
ode_observables{2, 1} = [30]; % Observable('tBid_', Bid(state=T))
ode_observables{2, 2} = [1];
ode_observables{2, 3} = 1;
ode_observables{3, 1} = [67 70]; % Observable('aSmac_', Smac(state=A))
ode_observables{3, 2} = [1 1];
ode_observables{3, 3} = 1;
ode_observables{4, 1} = [19 62 64]; % Observable('mSmac_', Smac(state=M))
ode_observables{4, 2} = [1 1 1];
ode_observables{4, 3} = 1;
ode_observables{5, 1} = [37]; % Observable('cPARP_', PARP(state=C))
ode_observables{5, 2} = [1];
ode_observables{5, 3} = 1;
end
