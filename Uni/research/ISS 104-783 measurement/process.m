close all; clear; clc; addpath(genpath(fullfile([pwd '\Uni' '\research' '\ISS 104-783 measurement'])))

%% initialization
c = 299792458;
MEASUREMENT = ["metal" "ceramics" "absorber"];

% ISS CPW parameters
ConductorWidth = 45e-6;         % measured
EpsilonR = 9.9;                 % nominal estimate (tune)
Height = 0.254e-3;              % nominal
SlotWidth = 30e-6;              % measured
Thickness = 3e-6;               % estimated (tune)

% models of calibration standards
ShortL = 2.4e-12;
OpenC = -9.3e-15;
MatchR = 47;                    % tune (nominal 50)
MatchL = -3.5e-12;              % tune (nominal -3.5e-12)

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

% Thru and Line{1..5} nominal delays from datasheet
% must be only for reference since propagation velocity is dispersive
% also the delays don't account for the skate mentioned above
ThruDelay = 1e-12;
Line1Delay = 3e-12;
Line2Delay = 7e-12;
Line3Delay = 14e-12;
Line4Delay = 27e-12;
Line5Delay = 40e-12;

% propagation velocity
LineVelocity = [
    ThruLength/ThruDelay
    Line1Length/Line1Delay
    Line2Length/Line2Delay
    Line3Length/Line3Delay
    Line4Length/Line4Delay
    Line5Length/Line5Delay
];
% eff. diel. const. for 1um thick gold plating using txline: 5.32 @ 30 GHz
LineEpsilon = (c./LineVelocity).^2;

% arbitrarily chosen measurement only to obtain frequency points
freq = cat(1, SXPParse(fullfile('metal', 'Open_DC_to_67_GHz.s2p')), SXPParse(fullfile('metal', 'Open_67_to_110_GHz.s2p')));

%% calibration and verification standards models
Z0 = 50;                        % measurement reference impedance

ShortImp = 1j*2*pi*freq*ShortL;
ShortRefl = (ShortImp - Z0)./(ShortImp + Z0);
OpenImp = 1./(1j*2*pi*freq*OpenC);
OpenRefl = (OpenImp - Z0)./(OpenImp + Z0);
MatchImp = MatchR + 1j*2*pi*freq*MatchL;
MatchRefl = (MatchImp - Z0)./(MatchImp + Z0);

% ports reflection coeficient
Port1Model = struct( ...
    'open', OpenRefl, ...
    'short', ShortRefl, ...
    'match', MatchRefl ...
);
Port2Model = Port1Model;

% CPW models of standards
ThruCPWModel = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', ThruLength ...
);
Line1CPWModel = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line1Length ...
);
Line2CPWModel = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line2Length ...
);
Line3CPWModel = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line3Length ...
);
Line4CPWModel = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line4Length ...
);
Line5CPWModel = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line5Length ...
);

% S-parameters of models
ThruModel = sparameters(ThruCPWModel, freq, MatchR);
Line1Model = sparameters(Line1CPWModel, freq, MatchR);
Line2Model = sparameters(Line2CPWModel, freq, MatchR);
Line3Model = sparameters(Line3CPWModel, freq, MatchR);
Line4Model = sparameters(Line4CPWModel, freq, MatchR);
Line5Model = sparameters(Line5CPWModel, freq, MatchR);

% imepdance of models (same for all lines)
[LineImpedance, EpsEff] = getZ0(ThruCPWModel, freq);

%% Phase error inspection on various pads for best pad material determination
% ax1 = matlab.graphics.axis.Axes.empty(3,0);
% ax2 = matlab.graphics.axis.Axes.empty(3,0);
% for m = 1:3
%     measurement = sprintf('%s', MEASUREMENT(m));
%     [~, P1Meas, P2Meas, ThruMeas, SLineMeas] = loadData(measurement);
% 
%     [Ea, Eb] = UOSM(P1Meas, P2Meas, ThruMeas, ...
%         Port1Model, Port2Model, squeeze(ThruModel.Parameters(2, 1, :)), freq);
% 
%     ThruCal = calibrate2PortSParDUT(Ea, Eb, ThruMeas);
%     Line1Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 1));
%     Line2Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 2));
%     Line3Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 3));
%     Line4Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 4));
%     Line5Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 5));
% 
%     ThruTransmPhaseErr = abs(unwrapPhase(angle(squeeze(ThruCal(2, 1, :)./ThruModel.Parameters(2, 1, :)))/pi*180));
%     Line1TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line1Cal(2, 1, :)./Line1Model.Parameters(2, 1, :)))/pi*180));
%     Line2TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line2Cal(2, 1, :)./Line2Model.Parameters(2, 1, :)))/pi*180));
%     Line3TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line3Cal(2, 1, :)./Line3Model.Parameters(2, 1, :)))/pi*180));
%     Line4TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line4Cal(2, 1, :)./Line4Model.Parameters(2, 1, :)))/pi*180));
%     Line5TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line5Cal(2, 1, :)./Line5Model.Parameters(2, 1, :)))/pi*180));
% 
%     figure(1)
%     ax1(m) = subplot(str2double(strcat(num2str(13), num2str(m))));
%     plot(freq/1e9, ThruTransmPhaseErr, 'Parent', ax1(m))
%     hold on
%     plot(freq/1e9, Line1TransmPhaseErr, 'Parent', ax1(m))
%     plot(freq/1e9, Line2TransmPhaseErr, 'Parent', ax1(m))
%     plot(freq/1e9, Line3TransmPhaseErr, 'Parent', ax1(m))
%     plot(freq/1e9, Line4TransmPhaseErr, 'Parent', ax1(m))
%     plot(freq/1e9, Line5TransmPhaseErr, 'Parent', ax1(m))
%     hold off
%     grid on
%     grid minor
%     xlabel('Frequency (GHz)')
%     ylabel('|\Delta(\angle S_{21})| (deg)')
%     linkaxes(ax1, 'xy')
%     legend(ax1(m), ...
%         strcat('Thru, mean: ', num2str(round(mean(ThruTransmPhaseErr), 2)), ', var: ', num2str(round(var(ThruTransmPhaseErr), 2))), ...
%         strcat('Line1, mean: ', num2str(round(mean(Line1TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line1TransmPhaseErr), 2))), ...
%         strcat('Line2, mean: ', num2str(round(mean(Line2TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line2TransmPhaseErr), 2))), ...
%         strcat('Line3, mean: ', num2str(round(mean(Line3TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line3TransmPhaseErr), 2))), ...
%         strcat('Line4, mean: ', num2str(round(mean(Line4TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line4TransmPhaseErr), 2))), ...
%         strcat('Line5, mean: ', num2str(round(mean(Line5TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line5TransmPhaseErr), 2))), ...
%         'Location', 'northwest' ...
%     )
%     title(['Measurement on ' measurement])
% 
%     figure(2)
%     ax2(m) = subplot(str2double(strcat(num2str(13), num2str(m))));
%     plot(freq/1e9, 20*log10(abs(squeeze(ThruCal(2, 1, :)))), 'Parent', ax2(m))
%     hold on
%     plot(freq/1e9, 20*log10(abs(squeeze(Line1Cal(2, 1, :)))), 'Parent', ax2(m))
%     plot(freq/1e9, 20*log10(abs(squeeze(Line2Cal(2, 1, :)))), 'Parent', ax2(m))
%     plot(freq/1e9, 20*log10(abs(squeeze(Line3Cal(2, 1, :)))), 'Parent', ax2(m))
%     plot(freq/1e9, 20*log10(abs(squeeze(Line4Cal(2, 1, :)))), 'Parent', ax2(m))
%     plot(freq/1e9, 20*log10(abs(squeeze(Line5Cal(2, 1, :)))), 'Parent', ax2(m))
%     yline(0)
%     hold off
%     grid on
%     grid minor
%     xlabel('Frequency (GHz)')
%     ylabel('|S_{21}| (dB)')
%     linkaxes(ax2, 'xy')
%     legend(ax2(m), 'Thru', 'Line1', 'Line2', 'Line3', 'Line4', 'Line5', 'Location', 'southwest')
%     title(['Measurement on ' measurement])
% end
% figure(1)
% sgtitle('Transmission phase error of calibrated lines')
% 
% figure(2)
% sgtitle('Transmission module of calibrated lines')

%% Calibration on metal: standards, error boxes
% m = 1;                          % choose metal for calibration
% measurement = sprintf('%s', MEASUREMENT(m));
% [~, P1Meas, P2Meas, ThruMeas, SLineMeas] = loadData(measurement);
% 
% [Ea, Eb] = UOSM(P1Meas, P2Meas, ThruMeas, ...
%     Port1Model, Port2Model, squeeze(ThruModel.Parameters(2, 1, :)), freq);
% 
% ThruCal = calibrate2PortSParDUT(Ea, Eb, ThruMeas);
% Line1Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 1));
% Line2Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 2));
% Line3Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 3));
% Line4Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 4));
% Line5Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 5));
% 
% figure(3)
% smithplot(freq, OpenRefl, 'LineWidth', 2, 'TitleTop', 'Calibration standards models')
% hold on
% smithplot(freq, ShortRefl, 'LineWidth', 2)
% smithplot(freq, MatchRefl, 'LineWidth', 2)
% legend('Open', 'Short', 'Match')
% 
% figure(4)
% showXPortParameters({freq}, {Ea, Eb}, 'Title', ['Error boxes ' measurement], ...
%    'LineWidth', [2 2], 'LegendText', {'Ea (S_{@@})', 'Eb (S_{@@})'})
% 
% figure(5)
% showXPortParameters({freq}, {ThruCal, ThruModel.Parameters}, 'Title', ['Calibrated Thru ' measurement], ...
%    'LineWidth', [2 2], 'LegendText', {'Meas (S_{@@})', 'Ideal (S_{@@})'})
% 
% figure(6)
% showXPortParameters({freq}, convert4DSparToCell(Line1Cal), 'Title', ['Calibrated Line 1 on ' measurement], 'Polar', true)
% 
% figure(7)
% showXPortParameters({freq}, convert4DSparToCell(Line2Cal), 'Title', ['Calibrated Line 2 on ' measurement], 'Polar', true)
% 
% figure(8)
% showXPortParameters({freq}, convert4DSparToCell(Line3Cal), 'Title', ['Calibrated Line 3 on ' measurement], 'Polar', true)
% 
% figure(9)
% showXPortParameters({freq}, convert4DSparToCell(Line4Cal), 'Title', ['Calibrated Line 4 on ' measurement], 'Polar', true)
% 
% figure(10)
% showXPortParameters({freq}, convert4DSparToCell(Line5Cal), 'Title', ['Calibrated Line 5 on ' measurement], 'Polar', true)

%% Results from the measurement on absorber for the tuning of match calibration standard
[~, P1Meas, P2Meas, ThruMeas, SLineMeas] = loadData('absorber');

[Ea, Eb] = UOSM(P1Meas, P2Meas, ThruMeas, ...
    Port1Model, Port2Model, squeeze(ThruModel.Parameters(2, 1, :)), freq);

ThruCal = calibrate2PortSParDUT(Ea, Eb, ThruMeas);
Line1Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 1));
Line2Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 2));
Line3Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 3));
Line4Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 4));
Line5Cal = calibrate2PortSParDUT(Ea, Eb, SLineMeas(:, :, :, 5));

% ThruTransmPhaseErr = abs(unwrapPhase(angle(squeeze(ThruCal(2, 1, :)./ThruModel.Parameters(2, 1, :)))/pi*180));
% Line1TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line1Cal(2, 1, :)./Line1Model.Parameters(2, 1, :)))/pi*180));
% Line2TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line2Cal(2, 1, :)./Line2Model.Parameters(2, 1, :)))/pi*180));
% Line3TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line3Cal(2, 1, :)./Line3Model.Parameters(2, 1, :)))/pi*180));
% Line4TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line4Cal(2, 1, :)./Line4Model.Parameters(2, 1, :)))/pi*180));
% Line5TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line5Cal(2, 1, :)./Line5Model.Parameters(2, 1, :)))/pi*180));
% 
% figure(11)
% plot(freq/1e9, ThruTransmPhaseErr)
% hold on
% plot(freq/1e9, Line1TransmPhaseErr)
% plot(freq/1e9, Line2TransmPhaseErr)
% plot(freq/1e9, Line3TransmPhaseErr)
% plot(freq/1e9, Line4TransmPhaseErr)
% plot(freq/1e9, Line5TransmPhaseErr)
% hold off
% grid on
% grid minor
% xlabel('Frequency (GHz)')
% ylabel('|\Delta(\angle S_{21})| (deg)')
% axis tight
% legend( ...
%     strcat('Thru, mean: ', num2str(round(mean(ThruTransmPhaseErr), 2)), ', var: ', num2str(round(var(ThruTransmPhaseErr), 2))), ...
%     strcat('Line1, mean: ', num2str(round(mean(Line1TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line1TransmPhaseErr), 2))), ...
%     strcat('Line2, mean: ', num2str(round(mean(Line2TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line2TransmPhaseErr), 2))), ...
%     strcat('Line3, mean: ', num2str(round(mean(Line3TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line3TransmPhaseErr), 2))), ...
%     strcat('Line4, mean: ', num2str(round(mean(Line4TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line4TransmPhaseErr), 2))), ...
%     strcat('Line5, mean: ', num2str(round(mean(Line5TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line5TransmPhaseErr), 2))), ...
%     'Location', 'northwest' ...
% )
% title('Transmission phase error of calibrated lines')

% figure(12)
% plot(freq/1e9, 20*log10(abs(squeeze(ThruCal(2, 1, :)))))
% hold on
% plot(freq/1e9, 20*log10(abs(squeeze(Line1Cal(2, 1, :)))))
% plot(freq/1e9, 20*log10(abs(squeeze(Line2Cal(2, 1, :)))))
% plot(freq/1e9, 20*log10(abs(squeeze(Line3Cal(2, 1, :)))))
% plot(freq/1e9, 20*log10(abs(squeeze(Line4Cal(2, 1, :)))))
% plot(freq/1e9, 20*log10(abs(squeeze(Line5Cal(2, 1, :)))))
% yline(0)
% hold off
% grid on
% grid minor
% xlabel('Frequency (GHz)')
% ylabel('|S_{21}| (dB)')
% axis tight
% legend('Thru', 'Line1', 'Line2', 'Line3', 'Line4', 'Line5', 'Location', 'southwest')
% title('Transmission module of calibrated lines')

figure(13)
plot(freq/1e9, 20*log10(abs(squeeze(ThruCal(1, 1, :)))))
hold on
plot(freq/1e9, 20*log10(abs(squeeze(Line1Cal(1, 1, :)))))
plot(freq/1e9, 20*log10(abs(squeeze(Line2Cal(1, 1, :)))))
plot(freq/1e9, 20*log10(abs(squeeze(Line3Cal(1, 1, :)))))
plot(freq/1e9, 20*log10(abs(squeeze(Line4Cal(1, 1, :)))))
plot(freq/1e9, 20*log10(abs(squeeze(Line5Cal(1, 1, :)))))
yline(0)
hold off
grid on
grid minor
xlabel('Frequency (GHz)')
ylabel('|S_{11}| (dB)')
axis tight
ylim([-70 0])
legend('Thru', 'Line1', 'Line2', 'Line3', 'Line4', 'Line5', 'Location', 'southeast')
title(['Reflection module of calibrated lines (MatchL = ' num2str(MatchL*1e12) ' pF)'])

% figure(14)
% plot(freq/1e9, LineImpedance)
% ylabel('Z (Ohm)')
% xlabel('Frequency (GHz)')
% title('CPW impedance according to txlineCPW')
% 
% figure(15)
% plot(freq/1e9, EpsEff)
% ylabel('\epsilon_{eff} (-)')
% xlabel('Frequency (GHz)')
% title('Effective permittivity of CPW according to txlineCPW')

%% Text outputs
LineVelocity = sprintf('%.2f, ', LineVelocity*1e-8);
LineVelocity = LineVelocity(1:end-2);
LineEpsilon = sprintf('%.2f, ', LineEpsilon);
LineEpsilon = LineEpsilon(1:end-2);

fprintf(['\n' ...
    '*******************************************\n' ...
    'STDOUT\n\n'])
fprintf(['Nominal parameters of lines based on the information in the ISS datasheet\n' ...
    '- the lines are ordered as {Thru, Line1, Line2, Line3, Line4, Line5}\n' ...
    '----\n'])
fprintf('LineVelocity (m/s):\t{%s} * 1e8\n', LineVelocity)
fprintf('LineEpsilon (-):\t{%s}\n', LineEpsilon)
