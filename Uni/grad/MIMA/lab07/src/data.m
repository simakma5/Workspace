close all; clear; clc

%% Task 2
S21 = [-14.98 -15.05 -15.22 -15.34 -15.57 -15.67 -15.76 -15.76 -14.75 -15.83];      % dB
S31 = [-3.347 -3.399 -3.524 -3.605 -3.717 -3.682 -3.705 -3.623 -3.553 -3.362];      % dB
P_DET = [-17.52 -17.57 -17.76 -17.78 -17.91 -17.87 -18.14 -18.31 -18.35 -18.80];    % dBm

S21 = 10.^(S21/20);                     % dB to linear
S31 = 10.^(S31/20);                     % dB to linear
P_DET = 10.^(P_DET/10)*1e-3;            % dBm to linear

P_LOAD = P_DET.*(S31./S21);             % calculation in linear scale
disp(round(10*log10(P_LOAD/1e-3),1))    % display in dBm

%% Task 4
c = 3e8;                                        % speed of light
D = 3.9e-2;                                     % antenna diagonal
F = [77.1 78 79 80 80.9];                       % GHz
Prx = [-25.59 -26.03 -24.73 -22.59 -23.41];     % dBm

d_F = 2*D^2*F(end)*1e9/c;
disp(['d_F = ' num2str(d_F*1e2) ' cm'])

d = 86.5e-2;                                    % cm, d > d_F
FSL = 20*log10(4*pi*d*F*1e9/c);
disp(round(FSL,1))

G = 25;                                         % dB, antenna gain
EIRP = Prx + FSL - G;
disp(round(EIRP,1))
