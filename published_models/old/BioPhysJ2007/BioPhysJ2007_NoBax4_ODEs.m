function [out] = BioPhysJ2007_NoBax4_ODEs(t, input, param)

Act = input(1);
Bcl2 = input(2);
AcBax = input(3);

k1 = param(1);
k2 = param(2);
k3 = param(3);
k4 = param(4);
k5 = param(5);
k6 = param(6);
k7 = param(7);
k8 = param(8);
k9 = param(9);
k10 = param(10);

Act_mito = param(11);
Bcl2_mito = param(12);
Bax_mito = param(13);

% Conservation Equations
ActBcl2 = Act_mito - Act;
AcBaxBcl2 = Bcl2_mito - ActBcl2 - Bcl2;
InBax = Bax_mito - AcBaxBcl2 - AcBax;

% Act
out(1,1) =  - k5 * Act * Bcl2 ...
            + k6 * ActBcl2 ...
            + k7 * AcBax ...
            - k8 * Act * AcBaxBcl2;
        
% AcBax
out(2,1) =  + k1 * Act * InBax ...
            - k2 * AcBax ...
            - k3 * AcBax * Bcl2 ...
            + k4 * AcBaxBcl2 ...
            - k7 * AcBax * ActBcl2 ...
            + k8 * Act * AcBaxBcl2;

% Bcl2
out(3,1) =  - k3 * AcBax ...
            - k5 * Act * Bcl2 ...
            + k6 * ActBcl2 ...
            + k4 * AcBaxBcl2;
end

