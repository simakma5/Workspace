clear
close all
clc

% task f
f0 = 10e9; % design frequency
subHeight = 0.508e-3; % height of the substrate
substrate = {'EpsilonR', 3.66, 'Height', subHeight, 'LossTangent', 0.0037, ...
   'Thickness', 17e-6, 'SigmaCond', 5.88e7}; % definition of the substrate
% width of microstrip line with impedance 50 Ohm
w = microstripW(substrate, f0, 50, 0.1, 0.1e-3, 10e-3, 1e-6);
nThick = 100; % number of thicknesses to check
% all line thicknesses, maximal thickness is subHeight/10 because of the limit
% in txlineMicrostrip function
lineThick = linspace(9e-6, subHeight/10, nThick);
% allocations
loss = zeros(nThick, 1);
for iThick = 1:nThick
   substrate{8} = lineThick(iThick); % change just the thickness
   % create single object representing the microstrip line
   microstrip = txlineMicrostrip(substrate{:}, 'LineLength', 1e-2, 'Width', w);
   % get S-parameters at frequency f0 with reference impedance 50 Ohm
   SObj = sparameters(microstrip, f0, 50);
   loss(iThick) = -20*log10(abs(SObj.Parameters(2, 1))); % computation of loss in dB
end

figure
plot(lineThick*1e6, loss)
xlabel('t (\mum)')
ylabel('Loss (dB/cm)')