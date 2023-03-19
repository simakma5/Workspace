clear
close all
clc

freqList = linspace(0.1, 20, 100)*1e9; % vector of frequencies
substrate = {'EpsilonR', 3.66, 'Height', 0.508e-3, 'LossTangent', 0.0037, ...
   'Thickness', 17e-6, 'SigmaCond', 5.88e7}; % definition of the substrate
w = 0.8e-3; % width from AWR TXLine
g = 0.2e-3; % gap from AWR TXLine
% create single object representing the CPWG line
CPWG = txlineCPW(substrate{:}, 'ConductorWidth', w, 'SlotWidth', g, ...
   'LineLength', 1e-2, 'ConductorBacked', true);
% get S-parameters of CPWG with reference impedance 50 Ohm
SObj = sparameters(CPWG, freqList, 50);
% get characteristic impedance
Z = getZ0(CPWG, freqList);
% plot results
figure
rfplot(SObj, 'mag')

figure
rfplot(SObj, 2, 1, 'angle')

figure
plot(freqList/1e9, Z)
xlabel('freq (GHz)')
ylabel('Z (\Omega)')