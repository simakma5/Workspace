close all; clear; clc; addpath(genpath(fullfile([pwd '\Uni' '\research' '\ISS 104-783 measurement'])))
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%% Initialization
% calibration Thru measured physical length
ThruLength = 200e-6;

% verification Line{1..5} nominal physical lengths from datasheet
Line1Length = 446e-6;
Line2Length = 896e-6;
Line3Length = 1796e-6;
Line4Length = 3496e-6;
Line5Length = 5246e-6;

% correction of lengths due to the vertical overtravel of probes during
% measurement resulting in a horizontal skate of probes
% - overtravel value recommended in the ISS datasheet: 75-125 um
% - skate value recommended in the ISS datasheet: 33 um
skate = 35e-6;
ThruLength =  ThruLength - 2*skate;
Line1Length = Line1Length - 2*skate;
Line2Length = Line2Length - 2*skate;
Line3Length = Line3Length - 2*skate;
Line4Length = Line4Length - 2*skate;
Line5Length = Line5Length - 2*skate;

% Line lengths minus Thru length
L = [ThruLength; Line1Length; Line2Length; Line3Length; Line4Length; Line5Length] - ThruLength;

%% Measurement data
% Thru measurement
[freq, ThruMeas] = SXPParse(fullfile('absorber', 'Thru_A1_DC_to_67_GHz.s2p'));

% Reflect measurement and parameters
[~, OpenMeas] = SXPParse(fullfile('absorber', 'Open_DC_to_67_GHz.s2p'));
[~, ShortMeas] = SXPParse(fullfile('absorber', 'Short_A1_DC_to_67_GHz.s2p'));
ReflectMeas = ShortMeas;
ReflectModel = -1;
ReflectOffset = -ThruLength/2;      % reference plane in the center of Thru

% Lines measurement
[~, Line1Meas] = SXPParse(fullfile('absorber', 'Line_A1_DC_to_67_GHz.s2p'));
[~, Line2Meas] = SXPParse(fullfile('absorber', 'Line_A2_DC_to_67_GHz.s2p'));
[~, Line3Meas] = SXPParse(fullfile('absorber', 'Line_A3_DC_to_67_GHz.s2p'));
[~, Line4Meas] = SXPParse(fullfile('absorber', 'Line_A4_DC_to_67_GHz.s2p'));
[~, Line5Meas] = SXPParse(fullfile('absorber', 'Line_A5_DC_to_67_GHz.s2p'));
LinesMeas = cat(4, ThruMeas, Line1Meas, Line2Meas, Line3Meas, Line4Meas, Line5Meas);

% Match measurement for calibration verifiaction
[~, MatchMeas] = SXPParse(fullfile('absorber', 'Match_A1_DC_to_67_GHz.s2p'));

%% TRL calibration
% error boxes for all Lines combinations
[allEa, allEb, propConst, EPD, lineComb, dL, reflectSol, epsEff] = ...
   mTRL(LinesMeas, L, ReflectMeas, ReflectModel, freq, 'm', ReflectOffset);

% Match and Open calibration with all combinations of error boxes
MatchAllCal = calibrate2PortSParDUT(allEa, allEb, MatchMeas);
OpenAllCal = calibrate2PortSParDUT(allEa, allEb, OpenMeas);

% weighted mean of all calibrated combinations
weights = sin(EPD).*(dL.').^2;
MatchCal = wMeanSPar(MatchAllCal, weights);
OpenCal = wMeanSPar(OpenAllCal, weights);

%% Figures
% legend for all combinations of Lines
lineCombLegend = sprintf('%d-%d\n', lineComb.');
lineCombLegend = strsplit(lineCombLegend, '\n');
lineCombLegend(end) = [];

% Reflect model for all line combinations (choose using the mouse wheel)
figure('Name', 'Calibrated Reflect Options')
inspectReflectSolutions(reflectSol, freq, lineComb, 'Calibrated Reflect Options')

% Error box A S-parameters
figure('Name', 'Error box A')
showXPortParameters({freq}, convert4DSparToCell(allEa), 'Title', 'Error box A')

% Error box B S-parameters
figure('Name', 'Error box B')
showXPortParameters({freq}, convert4DSparToCell(allEb), 'Title', 'Error box B')

% Propagation constant
% the longer the line difference, the more precise results we obtain
plotPropagationConst(freq, propConst, lineComb, dL, 'Propagation constant')

% Effective Phase Delay
figure('Name', 'Effective Phase Delay')
plot(freq*1e-9, EPD)
title('Effective Phase Delay')
grid on
xlabel('freq (GHz)')
ylabel('EPD')
legend(lineCombLegend)

% Effective Permittivity
figure('Name', 'Effective Permittivity')
plot(freq*1e-9, epsEff)
title('Effective Permittivity')
grid on
xlabel('freq (GHz)')
ylabel('\epsilon_{eff}')
legend(lineCombLegend)
ylim([0 10])

% Calibrated Match S-parameters for all Line combinations
figure('Name', 'Calibrated Match - all combinations')
showXPortParameters({freq}, convert4DSparToCell(MatchAllCal), ...
    'Title', 'Calibrated Match - all combinations')

% Calibrated Match S-parameters for the weighted average
figure('Name', 'Calibrated Match - weighted mean')
showXPortParameters({freq}, {MatchCal}, ...
    'Title', 'Calibrated Match - weighted mean')

% Calibrated Open S-parameters for all Line combinations
figure('Name', 'Calibrated Open - all combinations')
showXPortParameters({freq}, convert4DSparToCell(OpenAllCal), ...
    'Title', 'Calibrated Open - all combinations')

% Calibrated Open S-parameters for the weighted average
figure('Name', 'Calibrated Open - weighted mean')
showXPortParameters({freq}, {OpenCal}, ...
    'Title', 'Calibrated Open - weighted mean')
