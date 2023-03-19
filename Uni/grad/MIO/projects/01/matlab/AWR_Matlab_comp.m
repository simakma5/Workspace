clear
close all
clc

freqList = linspace(0.1, 20, 100)*1e9; % vector of frequencies
substrate = {'EpsilonR', 3.66, 'Height', 0.508e-3, 'LossTangent', 0.0037, ...
   'Thickness', 17e-6, 'SigmaCond', 5.88e7}; % definition of the substrate
w = 1.1e-3; % width from AWR TXLine
% create single object representing the microstrip line
microstrip = txlineMicrostrip(substrate{:}, 'Width', w, 'LineLength', 1e-2);
% get S-parameters of microstrip with reference impedance 50 Ohm
SObj = sparameters(microstrip, freqList, 50);
% write Touchstone file
rfwrite(SObj, 'microstrip_Matlab.s2p');
% plot results
figure
rfplot(SObj, 'mag')

figure
rfplot(SObj, 2, 1, 'angle')