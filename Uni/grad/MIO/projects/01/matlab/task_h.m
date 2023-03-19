clear
close all
clc

% task h
freqList = linspace(0.1, 20, 100)*1e9; % vector of frequencies
f0 = 10e9; % design frequency
substrate = {'EpsilonR', 3.66, 'Height', 0.508e-3, 'LossTangent', 0.0037, ...
   'Thickness', 17e-6, 'SigmaCond', 5.88e7}; % definition of the substrate
% width of microstrip line with 50 Ohm impedance
w = microstripW(substrate, f0, 50, 0.1, 0.1e-3, 10e-3, 1e-6);
% create single object representing the microstrip line
microstrip = txlineMicrostrip(substrate{:}, 'Width', w);
% get its char. impedance and effective permittivity at all frequencies
[Z, epsEff] = getZ0(microstrip, freqList);

figure
plot(freqList/1e9, Z)
xlabel('f (GHz)')
ylabel('Z (\Omega)')

figure
plot(freqList/1e9, epsEff)
xlabel('f (GHz)')
ylabel('\epsilon_{eff} (-)')