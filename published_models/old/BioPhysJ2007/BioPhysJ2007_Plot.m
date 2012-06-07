% Calculate amounts of each at mitochondria at start
Act_initial = 0;
Bax_initial = 2e-7;
Bcl2_initial = 1e-7;
Act_cytosol = 1e-7;
Bax_cytosol = 2e-7;
Bcl2_cytosol = 0;

F1 = 0;             % Bax
F2 = 1;             % Act
F3 = 0;             % Bcl2

ss_AcBax_norms = [];

tcs = figure;
for i = 1:length(F2);    
    Act_mito = Act_initial + (Act_cytosol * F2(i));
    Bcl2_mito = Bcl2_initial + (Bcl2_cytosol * F3);
    Bax_mito = Bax_initial + (Bax_cytosol * F1);

    [t x] = BioPhysJ2007_Run(Act_mito, Bcl2_mito, Bax_mito, 0);
    ss_AcBax_norm = (x(:,2)/Bax_mito);
    semilogx(t, ss_AcBax_norm);
    hold on;
    ss_AcBax_norms = [ss_AcBax_norms ss_AcBax_norm(end)]
end

figure;
plot(F2, ss_AcBax_norms);
title('Fig5');

