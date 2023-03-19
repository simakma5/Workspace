clear
close all
clc

freqList = linspace(0.1, 20, 100)*1e9; % vector of frequencies
substrate = {'EpsilonR', 3.66, 'DielectricThickness', 0.508e-3, 'LossTangent', 0.0037, ...
   'Thickness', 17e-6, 'SigmaCond', 5.88e7}; % definition of the substrate
% NOTE: DielectricThickness property defines just the half of the stack thickness
w = 0.5e-3; % width from AWR TXLine
% create single object representing the stripline
stripline = txlineStripline(substrate{:}, 'Width', w, 'LineLength', 1e-2);
% get S-parameters of stripline with reference impedance 50 Ohm
SObj = sparameters(stripline, freqList, 50);
% get characteristic impedance
Z = getZ0(stripline, freqList);
% plot results
figure
rfplot(SObj, 'mag')

figure
rfplot(SObj, 2, 1, 'angle')

figure
plot(freqList/1e9, Z)
xlabel('freq (GHz)')
ylabel('Z (\Omega)')