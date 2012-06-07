function [t x] = BioPhysJ2007_Run(Act_mito, Bcl2_mito, Bax_mito, UseBax4)

% Parameters
param(1) = 0.5e-6; % Units of M^-1 s^-1
param(2) = 0.1;    % Units of s^-1
param(3) = 2e-6;
param(4) = 0.001;
param(5) = 3e-6;
param(6) = 0.04;
param(7) = 2e-6;
param(8) = 0;
param(9) = 2e-6;
param(10) = 0;

% Initial values for calculating conservation equations
param(11) = Act_mito;
param(12) = Bcl2_mito;
param(13) = Bax_mito;

% Initial conditions for ODEs
initvals(1) = Act_mito;
initvals(2) = Bcl2_mito;
initvals(3) = 0;           % At start, no active Bax

% Time span for integration
tspan = [0 5000];

% Integrator settings
options=odeset('AbsTol', 1e-22, 'RelTol', 2.22045e-14);

% Run the simulation
[t x] = ode23s(@BioPhysJ2007_NoBax4_ODEs, tspan, initvals, options, param);

end

