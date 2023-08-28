close all
clear
clc
addpath(genpath(fileparts(which(mfilename))))

%% initialization
c = 299792458;
MEASUREMENT = ["metal" "ceramics" "absorber"];

% ISS CPW parameters
ConductorWidth = 60e-6;         % measured 45 um
EpsilonR = 9.9;                 % nominal
Height = 0.254e-3;              % nominal
SlotWidth = 33e-6;              % measured 30 um
Thickness = 3e-6;               % estimated (tune)

% models of calibration standards
ShortL = 2.4e-12;
OpenC = -9.3e-15;
MatchR = 50;
MatchL = -3.5e-12;

% calibration Thru parameters
ThruDelay = 1e-12;              % nominal delay
ThruLength = 200e-6;            % measured physical length

% verification Line{1..5} nominal parameters from datasheet
Line1Delay = 3e-12;
Line1Length = 446e-6;
Line2Delay = 7e-12;
Line2Length = 896e-6;
Line3Delay = 14e-12;
Line3Length = 1796e-6;
Line4Delay = 27e-12;
Line4Length = 3496e-6;
Line5Delay = 40e-12;
Line5Length = 5246e-6;

% length correction by probes' vertical overtravel (recommended 75-125 um)
skate = 35e-6;                  % probes' horizontal skate (rec. 33 um)
ThruLength =  ThruLength - 2*skate;
Line1Length = Line1Length - 2*skate;
Line2Length = Line2Length - 2*skate;
Line3Length = Line3Length - 2*skate;
Line4Length = Line4Length - 2*skate;
Line5Length = Line5Length - 2*skate;

% propagation velocity
LineVelocity = [
    ThruLength/ThruDelay
    Line1Length/Line1Delay
    Line2Length/Line2Delay
    Line3Length/Line3Delay
    Line4Length/Line4Delay
    Line5Length/Line5Delay
];
% eff. diel. const. for 1um thick gold plating using txline: 5.32304 @ 30 GHz
LineEpsilon = (c./LineVelocity).^2;

% arbitrarily chosen measurement only to obtain frequency points
freq = SXPParse(fullfile('metal', 'Open.s2p'));

%% models of standards
ShortImp = 1j*2*pi*freq*ShortL;
ShortRefl = (ShortImp - MatchR)./(ShortImp + MatchR);
OpenImp = 1./(1j*2*pi*freq*OpenC);
OpenRefl = (OpenImp - MatchR)./(OpenImp + MatchR);
MatchImp = 50 + 1j*2*pi*freq*MatchL;
MatchRefl = (MatchImp - MatchR)./(MatchImp + MatchR);

% reflection coeficient of ports
Port1Model = struct( ...
    'open', OpenRefl, ...
    'short', ShortRefl, ...
    'match', MatchRefl ...
);
Port2Model = Port1Model;

% CPW models of standards
ThruModel = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', ThruLength ...
);
Line1Model = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line1Length ...
);
Line2Model = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line2Length ...
);
Line3Model = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line3Length ...
);
Line4Model = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line4Length ...
);
Line5Model = txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line5Length ...
);

% S-parameters of models
ThruModelSparam = sparameters(ThruModel, freq, MatchR);
Line1ModelSparam = sparameters(Line1Model, freq, MatchR);
Line2ModelSparam = sparameters(Line2Model, freq, MatchR);
Line3ModelSparam = sparameters(Line3Model, freq, MatchR);
Line4ModelSparam = sparameters(Line4Model, freq, MatchR);
Line5ModelSparam = sparameters(Line5Model, freq, MatchR);

% imepdance of models (same for all lines)
LineImpedance = getZ0(ThruModel);

%% calibration on metal -- phase error inspection for best material determination
% ax1 = matlab.graphics.axis.Axes.empty(3,0);
% ax2 = matlab.graphics.axis.Axes.empty(3,0);
% for m = 1:3
%     measurement = sprintf('%s', MEASUREMENT(m));
%     [~, P1MeasMetal, P2MeasMetal, ThruMeasMetal, SLineMeasMetal] = loadData(measurement);
% 
%     [EaMetal, EbMetal] = UOSM(P1MeasMetal, P2MeasMetal, ThruMeasMetal, ...
%         Port1Model, Port2Model, squeeze(ThruModelSparam.Parameters(2, 1, :)), freq);
% 
%     ThruCalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, ThruMeasMetal);
%     Line1CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 1));
%     Line2CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 2));
%     Line3CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 3));
%     Line4CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 4));
%     Line5CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 5));
% 
%     ThruTransmPhaseErr = abs(unwrapPhase(angle(squeeze(ThruCalMetal(2, 1, :)./ThruModelSparam.Parameters(2, 1, :)))/pi*180));
%     Line1TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line1CalMetal(2, 1, :)./Line1ModelSparam.Parameters(2, 1, :)))/pi*180));
%     Line2TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line2CalMetal(2, 1, :)./Line2ModelSparam.Parameters(2, 1, :)))/pi*180));
%     Line3TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line3CalMetal(2, 1, :)./Line3ModelSparam.Parameters(2, 1, :)))/pi*180));
%     Line4TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line4CalMetal(2, 1, :)./Line4ModelSparam.Parameters(2, 1, :)))/pi*180));
%     Line5TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line5CalMetal(2, 1, :)./Line5ModelSparam.Parameters(2, 1, :)))/pi*180));
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
%     title(['Measurement on ' measurement])
%     grid on
%     grid minor
%     xlabel('Frequency [GHz]')
%     ylabel('|\Delta(\angle S_{21})| [deg]')
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
%     hold off
% 
%     figure(2)
%     ax2(m) = subplot(str2double(strcat(num2str(13), num2str(m))));
%     plot(freq/1e9, 20*log10(abs(squeeze(ThruCalMetal(2, 1, :)))), 'Parent', ax2(m))
%     hold on
%     plot(freq/1e9, 20*log10(abs(squeeze(Line1CalMetal(2, 1, :)))), 'Parent', ax2(m))
%     plot(freq/1e9, 20*log10(abs(squeeze(Line2CalMetal(2, 1, :)))), 'Parent', ax2(m))
%     plot(freq/1e9, 20*log10(abs(squeeze(Line3CalMetal(2, 1, :)))), 'Parent', ax2(m))
%     plot(freq/1e9, 20*log10(abs(squeeze(Line4CalMetal(2, 1, :)))), 'Parent', ax2(m))
%     plot(freq/1e9, 20*log10(abs(squeeze(Line5CalMetal(2, 1, :)))), 'Parent', ax2(m))
%     yline(0)
%     title(['Measurement on ' measurement])
%     grid on
%     grid minor
%     xlabel('Frequency [GHz]')
%     ylabel('|S_{21}| [dB]')
%     linkaxes(ax2, 'xy')
%     legend(ax2(m), 'Thru', 'Line1', 'Line2', 'Line3', 'Line4', 'Line5', 'Location', 'southwest')
%     hold off
% end
% figure(1)
% hold on
% sgtitle('Transmission phase error of calibrated lines')
% hold off
% 
% figure(2)
% hold on
% sgtitle('Transmission module of calibrated lines')
% hold off

%% calibration on metal -- standards, error boxes (unused, out-of-date)
% m = 1;
% measurement = sprintf('%s', MEASUREMENT(m));
% [~, P1MeasMetal, P2MeasMetal, ThruMeasMetal, SLineMeasMetal] = loadData(measurement);
% 
% [EaMetal, EbMetal] = UOSM(P1MeasMetal, P2MeasMetal, ThruMeasMetal, Port1Model, Port2Model, squeeze(ThruModelSparam(2, 1, :)), freq);
% 
% ThruCalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, ThruMeasMetal);
% Line1CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 1));
% Line2CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 2));
% Line3CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 3));
% Line4CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 4));
% Line5CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 5));
% 
% figure(3)
% smithplot(freq, OpenRefl, 'LineWidth', 2, 'TitleTop', 'Calibration standards models')
% hold on
% smithplot(freq, ShortRefl, 'LineWidth', 2)
% smithplot(freq, MatchRefl, 'LineWidth', 2)
% legend('Open', 'Short', 'Match')
% 
% figure(4)
% showXPortParameters({freq}, {EaMetal, EbMetal}, 'Title', ['Error boxes ' measurement], ...
%    'LineWidth', [2 2], 'LegendText', {'Ea (S_{@@})', 'Eb (S_{@@})'})
% 
% figure(5)
% showXPortParameters({freq}, {ThruCalMetal, ThruModelSparam}, 'Title', ['Calibrated Thru ' measurement], ...
%    'LineWidth', [2 2], 'LegendText', {'Meas (S_{@@})', 'Ideal (S_{@@})'})
% 
% figure(6)
% showXPortParameters({freq}, convert4DSparToCell(Line1CalMetal), 'Title', ['Calibrated Line 1 on ' measurement], 'Polar', true)
% 
% figure(7)
% showXPortParameters({freq}, convert4DSparToCell(Line2CalMetal), 'Title', ['Calibrated Line 2 on ' measurement], 'Polar', true)
% 
% figure(8)
% showXPortParameters({freq}, convert4DSparToCell(Line3CalMetal), 'Title', ['Calibrated Line 3 on ' measurement], 'Polar', true)
% 
% figure(9)
% showXPortParameters({freq}, convert4DSparToCell(Line4CalMetal), 'Title', ['Calibrated Line 4 on ' measurement], 'Polar', true)
% 
% figure(10)
% showXPortParameters({freq}, convert4DSparToCell(Line5CalMetal), 'Title', ['Calibrated Line 5 on ' measurement], 'Polar', true)

%% results from measurement on absorber for reference match tuning
[~, P1MeasMetal, P2MeasMetal, ThruMeasMetal, SLineMeasMetal] = loadData('absorber');

[EaMetal, EbMetal] = UOSM(P1MeasMetal, P2MeasMetal, ThruMeasMetal, ...
    Port1Model, Port2Model, squeeze(ThruModelSparam.Parameters(2, 1, :)), freq);

ThruCalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, ThruMeasMetal);
Line1CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 1));
Line2CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 2));
Line3CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 3));
Line4CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 4));
Line5CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 5));

ThruTransmPhaseErr = abs(unwrapPhase(angle(squeeze(ThruCalMetal(2, 1, :)./ThruModelSparam.Parameters(2, 1, :)))/pi*180));
Line1TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line1CalMetal(2, 1, :)./Line1ModelSparam.Parameters(2, 1, :)))/pi*180));
Line2TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line2CalMetal(2, 1, :)./Line2ModelSparam.Parameters(2, 1, :)))/pi*180));
Line3TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line3CalMetal(2, 1, :)./Line3ModelSparam.Parameters(2, 1, :)))/pi*180));
Line4TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line4CalMetal(2, 1, :)./Line4ModelSparam.Parameters(2, 1, :)))/pi*180));
Line5TransmPhaseErr = abs(unwrapPhase(angle(squeeze(Line5CalMetal(2, 1, :)./Line5ModelSparam.Parameters(2, 1, :)))/pi*180));

figure(11)
plot(freq/1e9, ThruTransmPhaseErr)
hold on
plot(freq/1e9, Line1TransmPhaseErr)
plot(freq/1e9, Line2TransmPhaseErr)
plot(freq/1e9, Line3TransmPhaseErr)
plot(freq/1e9, Line4TransmPhaseErr)
plot(freq/1e9, Line5TransmPhaseErr)
title('Transmission phase error of calibrated lines')
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|\Delta(\angle S_{21})| [deg]')
legend( ...
    strcat('Thru, mean: ', num2str(round(mean(ThruTransmPhaseErr), 2)), ', var: ', num2str(round(var(ThruTransmPhaseErr), 2))), ...
    strcat('Line1, mean: ', num2str(round(mean(Line1TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line1TransmPhaseErr), 2))), ...
    strcat('Line2, mean: ', num2str(round(mean(Line2TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line2TransmPhaseErr), 2))), ...
    strcat('Line3, mean: ', num2str(round(mean(Line3TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line3TransmPhaseErr), 2))), ...
    strcat('Line4, mean: ', num2str(round(mean(Line4TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line4TransmPhaseErr), 2))), ...
    strcat('Line5, mean: ', num2str(round(mean(Line5TransmPhaseErr), 2)), ', var: ', num2str(round(var(Line5TransmPhaseErr), 2))), ...
    'Location', 'northwest' ...
)
hold off

figure(12)
plot(freq/1e9, 20*log10(abs(squeeze(ThruCalMetal(2, 1, :)))))
hold on
plot(freq/1e9, 20*log10(abs(squeeze(Line1CalMetal(2, 1, :)))))
plot(freq/1e9, 20*log10(abs(squeeze(Line2CalMetal(2, 1, :)))))
plot(freq/1e9, 20*log10(abs(squeeze(Line3CalMetal(2, 1, :)))))
plot(freq/1e9, 20*log10(abs(squeeze(Line4CalMetal(2, 1, :)))))
plot(freq/1e9, 20*log10(abs(squeeze(Line5CalMetal(2, 1, :)))))
yline(0)
title('Transmission module of calibrated lines')
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|S_{21}| [dB]')
legend('Thru', 'Line1', 'Line2', 'Line3', 'Line4', 'Line5', 'Location', 'southwest')
hold off

%% Text outputs
format shortG
disp(' ')
disp('*******************************************')
disp('STDOUT')
disp(' ')

disp('Lines (CPW models of ISS Thru and Line{1..5}) parameters')
disp('----')
disp('LineVelocity [m/s]: 1e8 *')
disp(round(LineVelocity'./1e8, 2))
disp('LineEpsilon [-]:')
disp(round(LineEpsilon', 2))
disp(['LineImpedance [Ohm]: ' num2str(round(LineImpedance, 2))])