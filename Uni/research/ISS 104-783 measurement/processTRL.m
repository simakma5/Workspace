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
% - overtravel value recommended in the ISS datasheet: 75 to 125 um
% - skate value recommended in the ISS datasheet: 33 um
skate = 35e-6;
ThruLength =  ThruLength - 2*skate;
Line1Length = Line1Length - 2*skate;
Line2Length = Line2Length - 2*skate;
Line3Length = Line3Length - 2*skate;
Line4Length = Line4Length - 2*skate;
Line5Length = Line5Length - 2*skate;

% Line lengths minus Thru length
LinesLength = [ThruLength; Line1Length; Line2Length; Line3Length; Line4Length; Line5Length] - ThruLength;

%% Measurement data
% Thru measurement
folderName = 'absorber';
[freq67, ThruMeas67] = SXPParse(fullfile(folderName, 'Thru_A1_DC_to_67_GHz.s2p'));
[freq110, ThruMeas110] = SXPParse(fullfile(folderName, 'Thru_A2_67_to_110_GHz.s2p'));
freq = [freq67; freq110];

% Reflect measurement and parameters
[~, ShortMeas67] = SXPParse(fullfile(folderName, 'Short_A1_DC_to_67_GHz.s2p'));
[~, ShortMeas110] = SXPParse(fullfile(folderName, 'Short_A2_67_to_110_GHz.s2p'));

[~, OpenMeas67] = SXPParse(fullfile(folderName, 'Open_DC_to_67_GHz.s2p'));
[~, OpenMeas110] = SXPParse(fullfile(folderName, 'Open_67_to_110_GHz.s2p'));

ReflectModel = -1;              % using short as the reflect caliber
ReflectOffset = -ThruLength/2;  % reference plane in the center of Thru

% Lines measurement
[~, Line1Meas67] = SXPParse(fullfile(folderName, 'Line_A1_DC_to_67_GHz.s2p'));
[~, Line1Meas110] = SXPParse(fullfile(folderName, 'Line_A1_67_to_110_GHz.s2p'));

[~, Line2Meas67] = SXPParse(fullfile(folderName, 'Line_A2_DC_to_67_GHz.s2p'));
[~, Line2Meas110] = SXPParse(fullfile(folderName, 'Line_A2_67_to_110_GHz.s2p'));

[~, Line3Meas67] = SXPParse(fullfile(folderName, 'Line_A3_DC_to_67_GHz.s2p'));
[~, Line3Meas110] = SXPParse(fullfile(folderName, 'Line_A3_67_to_110_GHz.s2p'));

[~, Line4Meas67] = SXPParse(fullfile(folderName, 'Line_A4_DC_to_67_GHz.s2p'));
[~, Line4Meas110] = SXPParse(fullfile(folderName, 'Line_A4_67_to_110_GHz.s2p'));

[~, Line5Meas67] = SXPParse(fullfile(folderName, 'Line_A5_DC_to_67_GHz.s2p'));
[~, Line5Meas110] = SXPParse(fullfile(folderName, 'Line_A5_67_to_110_GHz.s2p'));

LineMeas67 = cat(4, ThruMeas67, Line1Meas67, Line2Meas67, Line3Meas67, Line4Meas67, Line5Meas67);
LineMeas110 = cat(4, ThruMeas110, Line1Meas110, Line2Meas110, Line3Meas110, Line4Meas110, Line5Meas110);

% Match measurement for calibration verifiaction
[~, MatchMeas67] = SXPParse(fullfile(folderName, 'Match_A1_DC_to_67_GHz.s2p'));
[~, MatchMeas110] = SXPParse(fullfile(folderName, 'Match_A2_67_to_110_GHz.s2p'));

%% TRL calibration
% error boxes for all Line combinations
% 
% The whole calibration process operates on measurements taken by an
% already pre-calibrated setup. This pre-calibration is performed at the
% tips of the coaxial cables connecting the VNA to the microwave probe
% station. However, the subsequent transmission lines leading to the probes
% are different in the measurement setups of the two frequency bands:
% # DC to 67 GHz: connectors, probes
% # 67 to 110 GHz: connectors, probes frequency adapters, WR-10 waveguides,
%   1mm coaxial cables
%
% Therefore the error boxes must be computed separately for both bands DC
% to 67 GHz and 67 to 110 GHz due to the fact that the two measurements are
% calibrated at completely different planes meaning that they calibrate
% different microwave circuits. Subsequently plotting the error boxes for
% the two bands connected makes little sense. Whereas the magnitude of
% resulting error-box S-parameters might connect relatively smoothly for a
% random reason, the phase will take on a completely different behaviour at
% 67 GHz and further.
[allEa67, allEb67, propConst67, EPD67, lineComb, dL, reflectSol67, epsEff67] = ...
   mTRL(LineMeas67, LinesLength, ShortMeas67, ReflectModel, freq67, 'm', ReflectOffset);
[allEa110, allEb110, propConst110, EPD110, ~, ~, reflectSol110, epsEff110] = ...
   mTRL(LineMeas110, LinesLength, ShortMeas110, ReflectModel, freq110, 'm', ReflectOffset);

% calibration for all Line combinations
MatchAllCal67 = calibrate2PortSParDUT(allEa67, allEb67, MatchMeas67);
MatchAllCal110 = calibrate2PortSParDUT(allEa110, allEb110, MatchMeas110);

OpenAllCal67 = calibrate2PortSParDUT(allEa67, allEb67, OpenMeas67);
OpenAllCal110 = calibrate2PortSParDUT(allEa110, allEb110, OpenMeas110);

ThruAllCal67 = calibrate2PortSParDUT(allEa67, allEb67, ThruMeas67);
ThruAllCal110 = calibrate2PortSParDUT(allEa110, allEb110, ThruMeas110);

% weighted mean of all calibration combinations
weights67 = sin(EPD67).*(dL.').^2;
weights110 = sin(EPD110).*(dL.').^2;

MatchCal = cat(3, wMeanSPar(MatchAllCal67, weights67), wMeanSPar(MatchAllCal110, weights110));
OpenCal = cat(3, wMeanSPar(OpenAllCal67, weights67), wMeanSPar(OpenAllCal110, weights110));
ThruCal = cat(3, wMeanSPar(ThruAllCal67, weights67), wMeanSPar(ThruAllCal110, weights110));

%% Plotting parameters
% legend for all Line combinations
combLegendText = sprintf('%d-%d (%3.1f mm)\n', [lineComb.'; dL.'*1e3]);
combLegendText = strsplit(combLegendText, '\n');
combLegendText(end) = [];

% indices of the combinations sorted by Line length
[~, ind] = sort(dL);

% combinations colors: greater Lines difference => hotter color
% (the greater the Lines difference, the more precise results we obtain)
allColors = mat2cell(jet(length(dL)), ones(length(dL), 1), 3);

%% Figures
% Reflect model for all Line combinations (choose using the mouse wheel)
inspectReflectSolutions([reflectSol67; reflectSol110], freq, lineComb, 'Calibrated Reflect Options')

% Error box A S-parameters
figure('Name', 'Error box A (DC to 67 GHz)')
showXPortParameters({freq67}, convert4DSparToCell(allEa67))

figure('Name', 'Error box A (67 to 110 GHz)')
showXPortParameters({freq110}, convert4DSparToCell(allEa110))

% Error box B S-parameters
figure('Name', 'Error box B (DC to 67 GHz)')
showXPortParameters({freq67}, convert4DSparToCell(allEb67))

figure('Name', 'Error box B (67 to 110 GHz)')
showXPortParameters({freq110}, convert4DSparToCell(allEb110))

% Propagation constant
plotPropagationConst(freq, [propConst67; propConst110], lineComb, dL, 'Propagation constant')

% Effective Phase Delay
figure('Name', 'Effective Phase Delay')
plot(freq/1e9, [EPD67; EPD110])
grid on
axis tight
xlabel('freq (GHz)')
ylabel('EPD')
legend(combLegendText, 'Location', 'best')
title('Effective Phase Delay')

% Effective Permittivity
figure('Name', 'Effective Permittivity')
hLinesEpsEff = plot(freq/1e9, [epsEff67(:, ind); epsEff110(:, ind)]);
[hLinesEpsEff.Color] = allColors{:};
grid on
axis tight
ylim([4 6])
xlabel('freq (GHz)')
ylabel('\epsilon_{eff} (-)')
legend(combLegendText(ind), 'Location', 'best')
title('Effective Permittivity')

% Calibrated Match S-parameters for all Line combinations
figure('Name', 'Calibrated Match - all Line combinations')
showXPortParameters({freq}, convert4DSparToCell(cat(3, MatchAllCal67, MatchAllCal110)))

% Calibrated Match S-parameters for the weighted average
figure('Name', 'Calibrated Match - weighted mean')
showXPortParameters({freq}, {MatchCal})

% Calibrated Open S-parameters for all Line combinations
figure('Name', 'Calibrated Open - all Line combinations')
showXPortParameters({freq}, convert4DSparToCell(cat(3, OpenAllCal67, OpenAllCal110)))

% Calibrated Open S-parameters for the weighted average
figure('Name', 'Calibrated Open - weighted mean')
showXPortParameters({freq}, {OpenCal})

% Calibrated Thru (the only known standard) S-parameters for all Line combinations
figure('Name', 'Calibrated Thru - all Line combinations')
showXPortParameters({freq}, convert4DSparToCell(cat(3, ThruAllCal67, ThruAllCal110)))

% Calibrated Thru (the only known standard) S-parameters for the weighted average
figure('Name', 'Calibrated Thru - weighted mean')
showXPortParameters({freq}, {ThruCal})

