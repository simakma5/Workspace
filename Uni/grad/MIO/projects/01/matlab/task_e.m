clear
close all
clc

% task e
f0 = 10e9;% design frequency
subHeight = 0.508e-3; % height of the substrate
substrate = {'EpsilonR', 3.66, 'Height', subHeight, 'LossTangent', 0.0037, ...
   'Thickness', 17e-6, 'SigmaCond', 5.88e7}; % definition of the substrate
ZMin = 2*4*pi*1e-7*subHeight*f0; % minimal impedance on single-mode propagation
% width of microstrip line with minimal impedance (maximal width)
wMax = microstripW(substrate, f0, ZMin, 0.1, 0.1e-3, 10e-3, 1e-6);
nWidths = 100; % number of widths to check
% all line widths, minimal 0.1 mm is fabrication limit
lineWidths = linspace(0.1e-3, wMax, nWidths);
% allocations
Z = zeros(nWidths, 1);
loss = zeros(nWidths, 1);
for iWidth = 1:nWidths
   % create single object representing the microstrip line
   microstrip = txlineMicrostrip(substrate{:}, 'LineLength', 1e-2, ...
      'Width', lineWidths(iWidth));
   % get its char. impedance at frequency f0
   Z(iWidth) = getZ0(microstrip, f0);
   % get S-parameters at frequency f0 with reference impedance as line's
   % impedance
   SObj = sparameters(microstrip, f0, Z(iWidth));
   loss(iWidth) = -20*log10(abs(SObj.Parameters(2, 1))); % computation od loss in dB
end

figure
plot(Z, loss)
xlabel('Z_v (\Omega)')
ylabel('Loss (dB/cm)')