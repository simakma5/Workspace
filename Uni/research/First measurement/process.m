clear
close all
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
    Line1Length/Line1Delay
    Line2Length/Line2Delay
    Line3Length/Line3Delay
    Line4Length/Line4Delay
    Line5Length/Line5Delay
];
% Eff. diel. const. for 1um thick gold plating using txline: 5.32304 @ 30 GHz
LineEpsilon = (c./LineVelocity).^2;

overtravel = 50e-6;             % estimate from microscope image
% effective parameters accounting for referance plane shift by overtravel
Line1Delay = Line1Delay - overtravel/LineVelocity(1);
Line1Length = Line1Length - overtravel;
Line2Delay = Line2Delay - overtravel/LineVelocity(2);
Line2Length = Line2Length - overtravel;
Line3Delay = Line3Delay - overtravel/LineVelocity(3);
Line3Length = Line3Length - overtravel;
Line4Delay = Line4Delay - overtravel/LineVelocity(4);
Line4Length = Line4Length - overtravel;
Line5Delay = Line5Delay - overtravel/LineVelocity(5);
Line5Length = Line5Length - overtravel;

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
), freq);
Line1Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line1Length ...
), freq);
Line2Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line2Length ...
), freq);
Line3Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line3Length ...
), freq);
Line4Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line4Length ...
), freq);
Line5Model = sparameters(txlineCPW( ...
    'ConductorWidth', ConductorWidth, ...
    'EpsilonR', EpsilonR, ...
    'Height', Height, ...
    'SlotWidth', SlotWidth, ...
    'Thickness', Thickness, ...
    'LineLength', Line5Length ...
), freq);

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
    
    figure(1)
    subplot(str2double(strcat(num2str(13), num2str(m))))
    plot(freq/1e9, unwrapPhase(angle(squeeze(ThruCalMetal(2, 1, :)./ThruModel.Parameters(2, 1, :)))/pi*180))
    hold on
    plot(freq/1e9, unwrapPhase(angle(squeeze(Line1CalMetal(2, 1, :)./Line1Model.Parameters(2, 1, :)))/pi*180))
    plot(freq/1e9, unwrapPhase(angle(squeeze(Line2CalMetal(2, 1, :)./Line2Model.Parameters(2, 1, :)))/pi*180))
    plot(freq/1e9, unwrapPhase(angle(squeeze(Line3CalMetal(2, 1, :)./Line3Model.Parameters(2, 1, :)))/pi*180))
    plot(freq/1e9, unwrapPhase(angle(squeeze(Line4CalMetal(2, 1, :)./Line4Model.Parameters(2, 1, :)))/pi*180))
    plot(freq/1e9, unwrapPhase(angle(squeeze(Line5CalMetal(2, 1, :)./Line5Model.Parameters(2, 1, :)))/pi*180))
    title(['Measurement on ' measurement])
    grid on
    grid minor
    xlabel('Frequency [GHz]')
    ylabel('\Delta(\angle S_{21}) [deg]')
    axis tight
    legend('Thru', 'Line1', 'Line2', 'Line3', 'Line4', 'Line5', 'Location', 'northwest')
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
