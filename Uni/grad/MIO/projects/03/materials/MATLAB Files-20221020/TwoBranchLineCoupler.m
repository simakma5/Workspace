clear
close all
clc

c0 = 299792458; % speed of light
% substrate definition
% NOTE: Thickness is just Height/10 because of model restrictions
substrate = {'EpsilonR', 3.66, 'Height', 0.127e-3, 'LossTangent', 0.0037, ...
   'Thickness', 12.7e-6, 'SigmaCond', 5.88e7}; 
f0 = 25e9; % design frequency

%% 50 Ohm lines
w50 = microstripW(substrate, f0, 50); % get width of 50 Ohm microtrip line
line1 = txlineMicrostrip(substrate{:}, 'Width', w50); % create microstrip object
[~, epsEff50] = getZ0(line1, f0); % get effective permittivity of the microstrip
lambda50 = c0/(f0*sqrt(epsEff50)); % wavelength of guided wave on misrostrip at f0
line1.LineLength = lambda50/4; % set length of microstrip line
line2 = clone(line1); % create secon 50 Ohm line (line2 = line1 does not work!)

%% 35 Ohm lines
w35 = microstripW(substrate, f0, 50/sqrt(2));
line3 = txlineMicrostrip(substrate{:}, 'Width', w35);
[~, epsEff35] = getZ0(line3, f0);
lambda35 = c0/(f0*sqrt(epsEff35));
line3.LineLength = lambda35/4;
line4 = clone(line3);

%% create a circuit
couplerObj = circuit();

couplerObj.add([1 2], line1); % connect line1 between nodes 1 and 2
couplerObj.add([3 4], line2);
couplerObj.add([1 3], line3);
couplerObj.add([2 4], line4);

% define node pairs as microwave ports, node 0 is a ground
couplerObj.setports([1 0], [2 0], [3 0], [4 0]);

%% analyze circuits
freq = linspace(0.1*f0, 2*f0); % vector of frequencies
SObj = sparameters(couplerObj, freq); % compute S-parameters

figure
rfplot(SObj, 1:4, 1) % same as rfplot(SObj, {[1 1]; [2 1]; [3 1]; [4 1]})

figure
rfplot(SObj, {[3 1]; [4 1]}, 'angle') % phase difference between S31 and S41 is 90 deg.