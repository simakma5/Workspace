close all
clear
clc

%% initialization
c = 299792458;
MEASUREMENT = ["metal" "ceramics" "absorber"];

% ISS CPW parameters
ConductorWidth = 45e-6;         % measured
EpsilonR = 9.9;                 % nominal
Height = 0.254e-3;              % nominal
SlotWidth = 30e-6;              % measured
Thickness = 1e-6;               % estimated (tune)

% models of calibration standards
ShortL = 2.4e-12;
OpenC = -9.3e-15;
MatchL = -3.5e-12;
MatchR = 50;

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

% effective parameters accounting for referance plane shift by overtravel
% ThruLength =  ThruLength - 95e-6;
% Line1Length = Line1Length - 95e-6;
% Line2Length = Line2Length - 105e-6;
% Line3Length = Line3Length - 110e-6;
% Line4Length = Line4Length - 130e-6;
% Line5Length = Line5Length - 140e-6;
ThruLength =  ThruLength - 100e-6;
Line1Length = Line1Length - 100e-6;
Line2Length = Line2Length - 100e-6;
Line3Length = Line3Length - 100e-6;
Line4Length = Line4Length - 100e-6;
Line5Length = Line5Length - 100e-6;

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
ThruModel = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', ThruLength ...
), freq, MatchR);
Line1Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line1Length ...
), freq, MatchR);
Line2Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line2Length ...
), freq, MatchR);
Line3Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line3Length ...
), freq, MatchR);
Line4Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line4Length ...
), freq, MatchR);
Line5Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line5Length ...
), freq, MatchR);

%% calibration on metal - phase error inspection
for m = 1:3
    measurement = sprintf('%s', MEASUREMENT(m));
    [~, P1MeasMetal, P2MeasMetal, ThruMeasMetal, SLineMeasMetal] = loadData(measurement);
    
    [EaMetal, EbMetal] = UOSM(P1MeasMetal, P2MeasMetal, ThruMeasMetal, ...
        Port1Model, Port2Model, squeeze(ThruModel.Parameters(2, 1, :)), freq);
    
    ThruCalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, ThruMeasMetal);
    Line1CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 1));
    Line2CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 2));
    Line3CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 3));
    Line4CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 4));
    Line5CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 5));
    
    ThruTransmPhaseErr = unwrapPhase(angle(squeeze(ThruCalMetal(2, 1, :)./ThruModel.Parameters(2, 1, :)))/pi*180);
    Line1TransmPhaseErr = unwrapPhase(angle(squeeze(Line1CalMetal(2, 1, :)./Line1Model.Parameters(2, 1, :)))/pi*180);
    Line2TransmPhaseErr = unwrapPhase(angle(squeeze(Line2CalMetal(2, 1, :)./Line2Model.Parameters(2, 1, :)))/pi*180);
    Line3TransmPhaseErr = unwrapPhase(angle(squeeze(Line3CalMetal(2, 1, :)./Line3Model.Parameters(2, 1, :)))/pi*180);
    Line4TransmPhaseErr = unwrapPhase(angle(squeeze(Line4CalMetal(2, 1, :)./Line4Model.Parameters(2, 1, :)))/pi*180);
    Line5TransmPhaseErr = unwrapPhase(angle(squeeze(Line5CalMetal(2, 1, :)./Line5Model.Parameters(2, 1, :)))/pi*180);

    figure(1)
    subplot(str2double(strcat(num2str(13), num2str(m))))
    plot(freq/1e9, ThruTransmPhaseErr)
    hold on
    plot(freq/1e9, Line1TransmPhaseErr)
    plot(freq/1e9, Line2TransmPhaseErr)
    plot(freq/1e9, Line3TransmPhaseErr)
    plot(freq/1e9, Line4TransmPhaseErr)
    plot(freq/1e9, Line5TransmPhaseErr)
    title(['Measurement on ' measurement])
    grid on
    grid minor
    xlabel('Frequency [GHz]')
    ylabel('\Delta(\angle S_{21}) [deg]')
    axis tight
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

    figure(2)
    subplot(str2double(strcat(num2str(13), num2str(m))))
    plot(freq/1e9, 20*log10(abs(squeeze(ThruCalMetal(2, 1, :)))))
    hold on
    plot(freq/1e9, 20*log10(abs(squeeze(Line1CalMetal(2, 1, :)))))
    plot(freq/1e9, 20*log10(abs(squeeze(Line2CalMetal(2, 1, :)))))
    plot(freq/1e9, 20*log10(abs(squeeze(Line3CalMetal(2, 1, :)))))
    plot(freq/1e9, 20*log10(abs(squeeze(Line4CalMetal(2, 1, :)))))
    plot(freq/1e9, 20*log10(abs(squeeze(Line5CalMetal(2, 1, :)))))
    yline(0)
    title(['Measurement on ' measurement])
    grid on
    grid minor
    xlabel('Frequency [GHz]')
    ylabel('|S_{21}| [dB]')
    axis tight
    legend('Thru', 'Line1', 'Line2', 'Line3', 'Line4', 'Line5', 'Location', 'southwest')
    hold off
end
figure(1)
hold on
sgtitle('Transmission phase error of calibrated lines')
hold off

figure(2)
hold on
sgtitle('Transmission module of calibrated lines')
hold off

disp(LineEpsilon')

%% calibration on metal - standards, error boxes, 
% m = 1;
% measurement = sprintf('%s', MEASUREMENT(m));
% [~, P1MeasMetal, P2MeasMetal, ThruMeasMetal, SLineMeasMetal] = loadData(measurement);
% 
% [EaMetal, EbMetal] = UOSM(P1MeasMetal, P2MeasMetal, ThruMeasMetal, Port1Model, Port2Model, squeeze(ThruModel(2, 1, :)), freq);
% 
% ThruCalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, ThruMeasMetal);
% Line1CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 1));
% Line2CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 2));
% Line3CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 3));
% Line4CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 4));
% Line5CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 5));
% 
% figure(2)
% smithplot(freq, OpenRefl, 'LineWidth', 2, 'TitleTop', 'Calibration standards models')
% hold on
% smithplot(freq, ShortRefl, 'LineWidth', 2)
% smithplot(freq, MatchRefl, 'LineWidth', 2)
% legend('Open', 'Short', 'Match')
% 
% figure(3)
% showXPortParameters({freq}, {EaMetal, EbMetal}, 'Title', ['Error boxes ' measurement], ...
%    'LineWidth', [2 2], 'LegendText', {'Ea (S_{@@})', 'Eb (S_{@@})'})
% 
% figure(4)
% showXPortParameters({freq}, {ThruCalMetal, ThruModel}, 'Title', ['Calibrated Thru ' measurement], ...
%    'LineWidth', [2 2], 'LegendText', {'Meas (S_{@@})', 'Ideal (S_{@@})'})
% 
% figure(5)
% showXPortParameters({freq}, convert4DSparToCell(Line1CalMetal), 'Title', ['Calibrated Line 1 on ' measurement], 'Polar', true)
% 
% figure(6)
% showXPortParameters({freq}, convert4DSparToCell(Line2CalMetal), 'Title', ['Calibrated Line 2 on ' measurement], 'Polar', true)
% 
% figure(7)
% showXPortParameters({freq}, convert4DSparToCell(Line3CalMetal), 'Title', ['Calibrated Line 3 on ' measurement], 'Polar', true)
% 
% figure(8)
% showXPortParameters({freq}, convert4DSparToCell(Line4CalMetal), 'Title', ['Calibrated Line 4 on ' measurement], 'Polar', true)
% 
% figure(9)
% showXPortParameters({freq}, convert4DSparToCell(Line5CalMetal), 'Title', ['Calibrated Line 5 on ' measurement], 'Polar', true)
