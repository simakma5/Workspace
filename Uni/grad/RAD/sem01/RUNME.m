clear
close all

clc
%% Surveillance radar parameters
f0 = 9.35e9; % carrier frequency
Pt = 2.5e3; % peak transmited power
tau = 2e-6; % width of pulse
PRI = 250e-6; % pulse repetition interval
T0 = 290; % noise temperature
G = 10^(30/10); % gain of an antenna
thetaAz = 2.4/180*pi; % 3 dB width of antenna pattern in azimuth
TSC = 1; % time of scanning 
NF = 10^(3/10); % noise figure of receiver
RCS = 0.1; % RCS of target

k = 1.38e-23; % Boltzmann constant

%% Achievable parameters
[dt, prf, pav, ep, ru] = pulseTrain(tau, PRI, Pt);

fprintf('Pulse repetition frequency: %.2f Hz\n', prf);
fprintf('Average transmitted power: %.2f W\n', pav);
fprintf('Pulse energy: %.2f J\n', ep);
fprintf('Unambiguous range: %.2f m\n', ru);

[delta_R] = rangeResolution(tau);

fprintf('Range resolution: %.2f m\n', delta_R);

R = 1e3;
Pr = receivedPower(Pt, f0, G, RCS, R);

fprintf('Received power: %.2f dBm\n', 10*log10(Pr/1e-3));

B = 1/tau;
Ni = k*T0*B;
fprintf('Noise power at antenna output: %.2f dBm\n', 10*log10(Ni/1e-3));

SNRi = Pr/Ni;
fprintf('SNR at antenna input: %.2f dB\n', 10*log10(SNRi)); % from single pulse

SNRiMin = 10^(10/10);
PrMin = Ni*NF*SNRiMin;
Rmax = maximalRange(Pt, G, f0, RCS, PrMin);

fprintf('Target maximal range: %.2f km\n', Rmax/1e3);

np = thetaAz*TSC*prf/(2*pi);

fprintf('Number of returnd pulses from single target: %.2f \n', np);
fprintf('Processing gain of coherent integration: %.2f dB\n', 10*log10(np));

PrMinI = PrMin/np;
RmaxI = maximalRange(Pt, G, f0, RCS, PrMinI);
fprintf('Target maximal range using integration: %.2f km\n', RmaxI/1e3);

%% Sweep some parameters

Rsw = linspace(1e3, 30e3).';
freq = [2e9, 5e9, 10e9];
Prsw = receivedPower(Pt, freq, G, RCS, Rsw);

figure
plot(Rsw/1e3, 10*log10(Prsw/1e-3))
legend('f=2e9', 'f=5e9', 'f=10e9')
grid on
ylabel('Received Power (dBm)')
xlabel('Target Distance (km)')


RCSsw = logspace(-2, 2, 5).';
PtSw = linspace(1, 5e3);
RmaxSw = maximalRange(PtSw, G, f0, RCSsw, PrMin);

figure
plot(PtSw, RmaxSw/1e3)
legendText = sprintf('RCS=%.2f,', RCSsw);
legendText = strsplit(legendText, ',');
legendText(end) = [];
legend(legendText{:})
grid on
ylabel('Maximal Target Range (km)')
xlabel('Transmitted Power (W)')


speedSw = linspace(0, 300).';
fd = dopplerFreq(freq, 0, speedSw);

figure
plot(speedSw, fd/1e3)
legend('f=2e9', 'f=5e9', 'f=10e9')
grid on
ylabel('Doppler Frequency (kHz)')
xlabel('Target Speed (m/s)')
