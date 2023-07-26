clear
close all
clc

%% initialization
MEASUREMENT = ["metal" "ceramics" "absorber"];

% models of calibration standards
LShort = 2.4e-12;
COpen = -9.3e-15;
LMatch = -3.5e-12;
RMatch = 50;

ThruLength = 200e-6;        % Thru length (physical, measured) 
ThruDelay = 1e-12;          % Thru delay (nominal)
Line1Delay = 3e-12;         % CPW Line 1 delay (nominal)
Line2Delay = 7e-12;         % CPW Line 2 delay (nominal)
Line3Delay = 14e-12;        % CPW Line 3 delay (nominal)
Line4Delay = 27e-12;        % CPW Line 4 delay (nominal)
Line5Delay = 40e-12;        % CPW Line 5 delay (nominal)

Z0 = 50;
% arbitrarily chosen measurement only to obtain frequency points
freq = SXPParse(fullfile('metal', 'Open.s2p'));
nFreq = length(freq);

%% models of standards
ShortImp = 1j*2*pi*freq*LShort;
ShortRefl = (ShortImp - Z0)./(ShortImp + Z0);
OpenImp = 1./(1j*2*pi*freq*COpen);
OpenRefl = (OpenImp - Z0)./(OpenImp + Z0);
MatchImp = 50 + 1j*2*pi*freq*LMatch;
MatchRefl = (MatchImp - Z0)./(MatchImp + Z0);

% reflection coeficient of ports
Port1Model = struct( ...
    'open', OpenRefl, ...
    'short', ShortRefl, ...
    'match', MatchRefl ...
);
Port2Model = Port1Model;

% model of Thru
cpwThru = txlineCPW( ...
    'ConductorWidth', 45e-6, ...
    'EpsilonR', 9.9, ...
    'Height', 0.254e-3, ...
    'SlotWidth', 30e-6, ...
    'Thickness', 1e-6, ...
    'LineLength', ThruLength ...
);
thruS21Model = squeeze(sparameters(cpwThru, freq).Parameters(2, 1, :));

% lossless models of Thru nad Line{1..5} based on their nominal delays
SIdealThru = transmLineSPar(shiftdim(1j*2*pi*freq*ThruDelay, -2), 1, Z0, Z0);
SIdealLine1 = transmLineSPar(shiftdim(1j*2*pi*freq*Line1Delay, -2), 1, Z0, Z0);
SIdealLine2 = transmLineSPar(shiftdim(1j*2*pi*freq*Line2Delay, -2), 1, Z0, Z0);
SIdealLine3 = transmLineSPar(shiftdim(1j*2*pi*freq*Line3Delay, -2), 1, Z0, Z0);
SIdealLine4 = transmLineSPar(shiftdim(1j*2*pi*freq*Line4Delay, -2), 1, Z0, Z0);
SIdealLine5 = transmLineSPar(shiftdim(1j*2*pi*freq*Line5Delay, -2), 1, Z0, Z0);

%% calibration on metal - phase error inspection
for m = 1:3
    measurement = sprintf('%s', MEASUREMENT(m));
    [~, P1MeasMetal, P2MeasMetal, ThruMeasMetal, SLineMeasMetal] = loadData(measurement);
    
    [EaMetal, EbMetal] = UOSM(P1MeasMetal, P2MeasMetal, ThruMeasMetal, ...
        Port1Model, Port2Model, squeeze(SIdealThru(2, 1, :)), freq);
    
    ThruCalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, ThruMeasMetal);
    Line1CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 1));
    Line2CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 2));
    Line3CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 3));
    Line4CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 4));
    Line5CalMetal = calibrate2PortSParDUT(EaMetal, EbMetal, SLineMeasMetal(:, :, :, 5));
    
    figure(1)
    subplot(str2double(strcat(num2str(13), num2str(m))))
    plot(freq/1e9, unwrapPhase(squeeze(angle(ThruCalMetal(2, 1, :)./SIdealThru(2, 1, :)))/pi*180))
    hold on
    plot(freq/1e9, unwrapPhase(squeeze(angle(Line1CalMetal(2, 1, :)./SIdealLine1(2, 1, :)))/pi*180))
    plot(freq/1e9, unwrapPhase(squeeze(angle(Line2CalMetal(2, 1, :)./SIdealLine2(2, 1, :)))/pi*180))
    plot(freq/1e9, unwrapPhase(squeeze(angle(Line3CalMetal(2, 1, :)./SIdealLine3(2, 1, :)))/pi*180))
    plot(freq/1e9, unwrapPhase(squeeze(angle(Line4CalMetal(2, 1, :)./SIdealLine4(2, 1, :)))/pi*180))
    plot(freq/1e9, unwrapPhase(squeeze(angle(Line5CalMetal(2, 1, :)./SIdealLine5(2, 1, :)))/pi*180))
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
    plot(freq/1e9, squeeze(20*log10(abs(ThruCalMetal(2, 1, :)))))
    hold on
    plot(freq/1e9, squeeze(20*log10(abs(Line1CalMetal(2, 1, :)))))
    plot(freq/1e9, squeeze(20*log10(abs(Line2CalMetal(2, 1, :)))))
    plot(freq/1e9, squeeze(20*log10(abs(Line3CalMetal(2, 1, :)))))
    plot(freq/1e9, squeeze(20*log10(abs(Line4CalMetal(2, 1, :)))))
    plot(freq/1e9, squeeze(20*log10(abs(Line5CalMetal(2, 1, :)))))
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
% [EaMetal, EbMetal] = UOSM(P1MeasMetal, P2MeasMetal, ThruMeasMetal, Port1Model, Port2Model, squeeze(SIdealThru(2, 1, :)), freq);
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
% showXPortParameters({freq}, {ThruCalMetal, SIdealThru}, 'Title', ['Calibrated Thru ' measurement], ...
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
