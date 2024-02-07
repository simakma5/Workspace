close all; clear; clc; addpath(genpath(fullfile(pwd, '\Uni', '\research', 'POTOM', '\ISS 104-783 measurement')))
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
c = 299792458;
epsilon0 = 1/(c^2*4*pi*1e-7);

%% Initialization

% ISS CPW parameters
W = 45e-6;                                                                  % measured conductor width
G = 30e-6;                                                                  % measured gap width
T = 3e-6;                                                                   % estimated conductor thickness (tune)
H = 0.254e-3;                                                               % nominal substrate height
epsilonR = 9.9;                                                             % nominal estimate of the substrate relative permittivity (tune)

% Calibration Thru measured physical length
ThruLength = 200e-6;

% Verification Line{1..5} nominal physical lengths from datasheet
Line1Length = 446e-6;
Line2Length = 896e-6;
Line3Length = 1796e-6;
Line4Length = 3496e-6;
Line5Length = 5246e-6;

% Correction of lengths due to the vertical overtravel of probes during
% measurement resulting in their horizontal skate over the contacts.
% The ISS datasheet recommends the overtravel of 75 to 125 um resulting
% in ~35 um of horizontal skate.
skate = 35e-6;
ThruLength =  ThruLength - 2*skate;
Line1Length = Line1Length - 2*skate;
Line2Length = Line2Length - 2*skate;
Line3Length = Line3Length - 2*skate;
Line4Length = Line4Length - 2*skate;
Line5Length = Line5Length - 2*skate;

% Lines length minus Thru length
LinesLength = [ThruLength; Line1Length; Line2Length; Line3Length; Line4Length; Line5Length] - ThruLength;

%% Measurement data
% Thru measurement
folderName = 'absorber';
[f67, ThruMeas67] = SXPParse(fullfile(folderName, 'Thru_A1_DC_to_67_GHz.s2p'));
[f110, ThruMeas110] = SXPParse(fullfile(folderName, 'Thru_A2_67_to_110_GHz.s2p'));
f = [f67; f110];
omega = 2*pi*f;

% Reflect measurement and parameters
[~, ShortMeas67] = SXPParse(fullfile(folderName, 'Short_A1_DC_to_67_GHz.s2p'));
[~, ShortMeas110] = SXPParse(fullfile(folderName, 'Short_A2_67_to_110_GHz.s2p'));

[~, OpenMeas67] = SXPParse(fullfile(folderName, 'Open_DC_to_67_GHz.s2p'));
[~, OpenMeas110] = SXPParse(fullfile(folderName, 'Open_67_to_110_GHz.s2p'));

ReflectModel = -1;                                                          % using short as the reflect caliber
ReflectOffset = -ThruLength/2;                                              % reference plane in the center of Thru

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
% The whole calibration process operates on measurements taken by an
% already pre-calibrated setup. This pre-calibration is performed at the
% tips of the coaxial cables connecting the VNA to the microwave probe
% station. However, the subsequent transmission lines leading to the probes
% are different in the measurement setups of the two frequency bands:
% * DC to 67 GHz: connectors, probes
% * 67 to 110 GHz: connectors, probes, frequency adapters, WR-10
%   waveguides, 1mm coaxial cables
%
% Therefore the error boxes must be computed separately for both bands DC
% to 67 GHz and 67 to 110 GHz due to the fact that the two measurements are
% calibrated at completely different planes meaning that they calibrate
% different microwave circuits. Subsequently plotting the error boxes for
% the two bands connected makes little sense. Whereas the magnitude of
% resulting error-box S-parameters might connect relatively smoothly for a
% random reason, the phase will take on a completely different behaviour at
% 67 GHz and further.

% Error boxes for all Line combinations
[allEa67, allEb67, propConst67, EPD67, lineComb, dL, reflectSol67, epsEff67] = ...
   mTRL(LineMeas67, LinesLength, ShortMeas67, ReflectModel, f67, 'm', ReflectOffset);
[allEa110, allEb110, propConst110, EPD110, ~, ~, reflectSol110, epsEff110] = ...
   mTRL(LineMeas110, LinesLength, ShortMeas110, ReflectModel, f110, 'm', ReflectOffset);
propConst = [propConst67; propConst110];

% Calibration for all Line combinations
MatchAllCal67 = calibrate2PortSParDUT(allEa67, allEb67, MatchMeas67);
MatchAllCal110 = calibrate2PortSParDUT(allEa110, allEb110, MatchMeas110);

OpenAllCal67 = calibrate2PortSParDUT(allEa67, allEb67, OpenMeas67);
OpenAllCal110 = calibrate2PortSParDUT(allEa110, allEb110, OpenMeas110);

ThruAllCal67 = calibrate2PortSParDUT(allEa67, allEb67, ThruMeas67);
ThruAllCal110 = calibrate2PortSParDUT(allEa110, allEb110, ThruMeas110);

% Weighted mean of all calibration combinations
weights67 = sin(EPD67).*(dL.').^2;
weights110 = sin(EPD110).*(dL.').^2;

MatchCal = cat(3, wMeanSPar(MatchAllCal67, weights67), wMeanSPar(MatchAllCal110, weights110));
OpenCal = cat(3, wMeanSPar(OpenAllCal67, weights67), wMeanSPar(OpenAllCal110, weights110));
ThruCal = cat(3, wMeanSPar(ThruAllCal67, weights67), wMeanSPar(ThruAllCal110, weights110));

%% In search of impedance (text)
% From [Marks and Williams, Characteristic Impedance Determination Using
% Propagation Constant Measurement] (paraphrased):
%
% Assumptions:
% * low substrate loss,
% * weak transverse current in the conductors.
% These are is typically true except at very high frequencies.
%
% This means $G$ is negligible and thus the equation
% $\gamma/Z_0 = i\omega C + G,$                                         (1)
% can be used for the determination of the phase of $Z_0$ from a
% measurement of the phase of $\gamma$, the
% propagation constant readily given by the TRL calibration. The allows
% for the transformation to a real impedance.
% If, in addition, $C$ is known then Equation~1 allows for the
% determination of the magnitude of $Z_0$ from $\gamma$ as well.
%
% Note: A similar way of characteristic impedance determination can be
% developed using the values of $L$ and $R$. However, this approach can not
% be recommended over the one described above due to the intrinsic problems
% caused by the strong dependence of $L$ and $R$ on the conductivity and
% frequency due to the varying current inside the metal.
%
%
% From [Lipa, Steer, Morris and Franzon, Comparison of Methods for
% Determining the Capacitance of Planar Transmission Lines with Application
% to Multichip Module Characterization] (paraphrased):
% 
% Method B (High Frequency Extrapolation):
% Assumed frequency dependence of the per unit components of the line:
% * $C$ and $G/\omega$ are approximately frequency independent,
% * $R$ has a DC component $R_{DC}$ and a skin effect resistance $R_S$
%   which has an assumed approximately $\omega^0.5$ dependence when the
%   skin effect is fully established,
% * $L = L_0 + L_int$, where $L_int$ is the internal conductor inductance
%   of the line due to current internal to the conductors and is
%   asymptotically zero at high frequencies (fully establised skin effect),
% * $R << \omega L_0$ at high frequencies.
% Result:
% $C \approx c^2*C_0/\omega^2 Re(\gamma^2),$                            (2)
% where $C_0$ is the per unit capacitance of the line in the absence of
% dielectric.
% Personal note: The approach described above assumes we know $C_0$.
% However, this knowledge can be utilized to make a direct estimate of the
% characteristic impedance without estimating the value of $C$. This
% estimate is then given as
% $Z_C = -i\gamma/(\omega\epsilon_{r,eff}C_0}.$                         (3)
% This, on the other hand, assumes we know $\epsilon_{r,eff}$. This
% quantity already conveys the frequency dependence of $C$.
%
% Method C (Low Frequency Extrapolation - $R_{DC}$ method):
% Assumptions:
% * negligible dielectric loss $LG$,
% * negligible skin effect resistance $R_S$.
% These simplifications are reasonable at low frequencies.
% Result:
% $C = \lim_{f \to f_0} Im(\gamma^2)/(R_{DC}*\omega).$                  (4)

% Indices of the Line combinations sorted by Line length difference
[~, ind] = sort(dL);
% Select the greatest length difference measurement for best accuracy
Gamma = propConst(:,ind(end));

%% Method B for capacitance per unit length determination
% Conformal maps parameters
k1 = W/(W+2*G);
k1Prime = sqrt(1-k1^2);
k2 = sinh(pi*W/(4*H))/sinh(pi*(W+2*G)/(4*H));
k2Prime = sqrt(1-k2^2);
% % Elliptic integral ratio approximation (error < 8e-6)
% if 0 <= k1 && k1 <= 1/sqrt(2)
%     EIRatio = pi/log(2*(1+sqrt(k1Prime))/(1-sqrt(k1Prime)));
% elseif 1/sqrt(2) <= k1 && k1 <= 1
%     EIRatio = log(2*(1+sqrt(k1))/(1-sqrt(k1)))/pi;
% else
%     error("The conformal-map variable k1 is out of the interval [0,1]!")
% end
% % Approximation verification: FAIL - error > 8e-6! ... discuss
% assert(EIRatio-ellipke(k1)/ellipke(k1Prime));
EIRatio1 = ellipke(k1)/ellipke(k1Prime);
EIRatio2 = ellipke(k2)/ellipke(k2Prime);
epsilonEff = 1+(epsilonR-1)*EIRatio2/(2*EIRatio1);

% DC capacitance per unit length (B2M17MIOA lectures from prof. Hoffmann)
C_t = epsilon0*(epsilonEff*4*EIRatio1+(2.3*1.65^(-18*T/G)+2)*T/G);
C_1t = epsilon0*(4*EIRatio1+(2.3*1.65^(-18*T/G)+2)*T/G);
C_0 = C_t + C_1t;                                                           % I guess that's how you get the final value?

% HF approximation of capacitance per unit length
C_HF = c^2*C_0*real(Gamma.^2)./(omega.^2);

% Characteristic impedance
Z_HF = -1i*Gamma./(omega.*C_HF);

%% Direct computation of impedance from $C_0$
Z_Direct = -1i*Gamma./(omega*epsilonEff*C_0);

%% LF extrapolation of capacitance per unit length ($R_{DC}$ method)
R_DC = 900e-3;                                                                 % Not exact so far
C_LF = imag(Gamma.^2)./(R_DC*omega);

% Characteristic impedance
Z_LF = -1i*Gamma./(omega.*C_LF);

%% Figures
% Legend for all Line combinations
combLegendText = sprintf('%d-%d (%3.1f mm)\n', [lineComb.'; dL.'*1e3]);
combLegendText = strsplit(combLegendText, '\n');
combLegendText(end) = [];

% Combinations colors: Lines difference ~ color temperature
% (the greater the Lines difference, the more precise results we obtain)
allColors = mat2cell(jet(length(dL)), ones(length(dL), 1), 3);

% Impedances concatenation for easier plotting
Z = [Z_Direct, Z_HF, Z_LF];

% % Reflect model for all Line combinations (choose using the mouse wheel)
% inspectReflectSolutions([reflectSol67; reflectSol110], f, lineComb, 'Calibrated Reflect Options')
% 
% % Error box A S-parameters
% figure('Name', 'Error box A (DC to 67 GHz)')
% showXPortParameters({freq67}, convert4DSparToCell(allEa67))
% 
% figure('Name', 'Error box A (67 to 110 GHz)')
% showXPortParameters({freq110}, convert4DSparToCell(allEa110))
% 
% % Error box B S-parameters
% figure('Name', 'Error box B (DC to 67 GHz)')
% showXPortParameters({freq67}, convert4DSparToCell(allEb67))
% 
% figure('Name', 'Error box B (67 to 110 GHz)')
% showXPortParameters({freq110}, convert4DSparToCell(allEb110))
% 
% % Propagation constant
% plotPropagationConst(f, propConst, lineComb, dL, 'Propagation constant');
% 
% % Effective Phase Delay
% figure('Name', 'Effective Phase Delay')
% plot(f/1e9, [EPD67; EPD110])
% grid on
% axis tight
% xlabel('Frequency (GHz)')
% ylabel('EPD')
% legend(combLegendText, 'Location', 'best')
% title('Effective Phase Delay')
% 
% % Effective Permittivity
% figure('Name', 'Effective Permittivity')
% hLinesEpsEff = plot(f/1e9, [epsEff67(:, ind); epsEff110(:, ind)]);
% [hLinesEpsEff.Color] = allColors{:};
% grid on
% axis tight
% ylim([4 6])
% xlabel('Frequency (GHz)')
% ylabel('\epsilon_{eff} (-)')
% legend(combLegendText(ind), 'Location', 'best')
% title('Effective Permittivity')
% 
% % Calibrated Match S-parameters for all Line combinations
% figure('Name', 'Calibrated Match - all Line combinations')
% showXPortParameters({f}, convert4DSparToCell(cat(3, MatchAllCal67, MatchAllCal110)))
% 
% % Calibrated Match S-parameters for the weighted average
% figure('Name', 'Calibrated Match - weighted mean')
% showXPortParameters({f}, {MatchCal})
% 
% % Calibrated Open S-parameters for all Line combinations
% figure('Name', 'Calibrated Open - all Line combinations')
% showXPortParameters({f}, convert4DSparToCell(cat(3, OpenAllCal67, OpenAllCal110)))
% 
% % Calibrated Open S-parameters for the weighted average
% figure('Name', 'Calibrated Open - weighted mean')
% showXPortParameters({f}, {OpenCal})
% 
% % Calibrated Thru (the only known standard) S-parameters for all Line combinations
% figure('Name', 'Calibrated Thru - all Line combinations')
% showXPortParameters({f}, convert4DSparToCell(cat(3, ThruAllCal67, ThruAllCal110)))
% 
% % Calibrated Thru (the only known standard) S-parameters for the weighted average
% figure('Name', 'Calibrated Thru - weighted mean')
% showXPortParameters({f}, {ThruCal})

% Characteristic impedance (real part)
figure('WindowState', 'maximized', 'Name', 'Characteristic impedance');
tiles = tiledlayout(2, 1);
nexttile;
plot(f/1e9, real(Z_HF))
hold on
plot(f/1e9, real(Z_Direct))
plot(f/1e9, real(Z_LF))
yline(50, '--', '50')
hold off
grid on
xlim([0 110])
ylim([-100 100])
% xlabel('Frequency [GHz]')
ylabel('R_0 [\Omega]')

nexttile;
plot(f/1e9, imag(Z_HF))
hold on
plot(f/1e9, imag(Z_Direct))
plot(f/1e9, imag(Z_LF))
hold off
grid on
xlim([0 110])
ylim([-10 10])
xlabel('Frequency [GHz]')
ylabel('X_0 [\Omega]')
legend( ...
    'High-frequency capacitance approximation', ...
    'Direct computation approximating capacitance as C_0', ...
    'Low-frequency capacitance extrapolation (R_{DC} method)', ...
    'Location', 'northoutside' ...
)
title(tiles, 'Characteristic impedance');
