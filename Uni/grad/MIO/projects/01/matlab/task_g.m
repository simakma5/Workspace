clear
close all
clc

% task e
f0 = 10e9;% design frequency
subHeight = 0.508e-3; % height of the substrate
substrate = {'EpsilonR', 3.66, 'Height', subHeight, 'LossTangent', 0.0037, ...
   'Thickness', 17e-6, 'SigmaCond', 5.88e7}; % definition of the substrate
ZMin = 2*4*pi*1e-7*subHeight*f0; % minimal impedance on single-mode propagation
wMin = 0.1e-3; % minimal width of the microstrip
% maximal impedance of the microstrip (minimal width0
ZMax = getZ0(txlineMicrostrip(substrate{:}, 'Width', wMin), f0);
% width of microstrip line with minimal impedance (maximal width)
wMax = microstripW(substrate, f0, ZMin, 0.1, 0.1e-3, 10e-3, 1e-6);
nThick = 100; % number of thicknesses to check
% all line thicknesses, maximal thickness is subHeight/10 because of the limit
% in txlineMicrostrip function
lineThick = linspace(9e-6, subHeight/10, nThick);
% allocations
ZMaxW = zeros(nThick, 1);
ZMinW = zeros(nThick, 1);
for iThick = 1:nThick
   substrate{8} = lineThick(iThick); % change just the thickness
   % create single object representing the microstrip line
   microstripMaxW = txlineMicrostrip(substrate{:}, 'Width', wMax);
   microstripMinW = txlineMicrostrip(substrate{:}, 'Width', wMin);
   % get its char. impedance and effectiove permittivity at frequency f0
   ZMaxW(iThick) = getZ0(microstripMaxW, f0);
   ZMinW(iThick) = getZ0(microstripMinW, f0);
end

figure
plot(lineThick*1e6, ZMaxW)
hold on
plot(lineThick*1e6, ZMinW)
xlabel('t (\mum)')
ylabel('Z_v (\Omega)')