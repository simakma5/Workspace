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
