%% Signal Detection in AWGN 

%% Single Sample Detection Using Coherent Receiver
% 100000-trial Monte-Carlo simulation

% fix the random number generator
rstream = RandStream.create('mt19937ar','seed',2009);

Ntrial = 1e5;             % number of Monte-Carlo trials
snrdb = 3;                % SNR in dB
snr = db2pow(snrdb);      % SNR in linear scale
spower = 1;               % signal power is 1
npower = spower/snr;           % noise power
namp = sqrt(npower/2);    % noise amplitude in each channel
s = ones(1,Ntrial);       % simple signal realization as constant value 
n = namp*(randn(rstream,1,Ntrial)+1i*randn(rstream,1,Ntrial));  % noise complex, white and Gaussian distributed
%%
%
% If the received signal contains the target, it is given by
x = s + n;

%%
% The matched filter in this case is trivial, since the signal itself is
% a unit sample.

mf = 1;

%% 
% In this case, the matched filter gain is 1, therefore, there is no SNR
% gain.  
%
% Now we do the detection and examine the performance of the detector. For
% a coherent receiver, the received signal after the matched filter is
% given by

y = mf'*x;  % apply the matched filter

%%
% The sufficient statistic, i.e., the value used to compare to the
% detection threshold, for a coherent detector is the real part of the
% received signal after the matched filter, i.e.,

z = real(y);

%%
%
% The required SNR threshold given a complex, white Gaussian noise for the
% NP detector for required Pfa can be calculated using the npwgnthresh function as follows:

Pfa = 1e-3;
% snrthreshold = db2pow(npwgnthresh1(Pfa, 1,'coherent'));
npulses = 1;
% snrthreshold = db2pow(mag2db(abs(erfcinv(2*Pfa)*sqrt(npulses))));
% SNR threshold for coherent case
snrthreshold = (abs(erfcinv(2*Pfa)*sqrt(npulses)))^2;

%%
% Note that this threshold, although also in the form of an SNR value, is
% different to the SNR of the received signal. The threshold SNR is a
% calculated value based on the desired detection performance, in this case
% the _Pfa_; while the received signal SNR is the physical characteristic
% of the signal determined by the propagation environment, the waveform,
% the transmit power, etc.
%
% The true threshold _T_ can then be derived from this SNR threshold as
%
% $$T=\sqrt{NM}\cdot\sqrt{\rm SNR}.$$
%

mfgain = mf'*mf;
% To match the equation in the text above
% npower - N
% mfgain - M
% snrthreshold - SNR
threshold = sqrt(npower*mfgain*snrthreshold);

%%
% The detection is performed by comparing the signal to the threshold.
% Since the original signal, _s_, is presented in the received signal, a
% successful detection occurs when the received signal passes the
% threshold, i.e. _z>T_. The capability of the detector to detect a target
% is often measured by the _Pd_. In a Monte-Carlo simulation, _Pd_ can be
% calculated as the ratio between the number of times the signal passes the
% threshold and the number of total trials.

Pd = sum(z>threshold)/Ntrial

%%
% On the other hand, a false alarm occurs when the detection shows that
% there is a target but there actually isn't one, i.e., the received signal
% passes the threshold when there is only noise present. The error
% probability of the detector to detect a target when there isn't one is
% given by _Pfa_. 

x = n;
y = mf'*x;
z = real(y);
Pfa = sum(z>threshold)/Ntrial

%%
% which meets our requirement.
%
% To see the relation among SNR, _Pd_ and _Pfa_ in a graph, we can plot the
% theoretical ROC curve using the rocsnr function for a SNR value
rocsnr1(snrdb,'SignalType','NonfluctuatingCoherent','MinPfa',1e-6);
SNR=snrdb;
MinPfa = 1e-6;
MaxPfa = 1;
NumPoints =101;
NumPulses=1;
d = db2pow(SNR);
MinPfa = log10(MinPfa);  % convert to log scale for even space sampling
MaxPfa = log10(MaxPfa);
Pfa = logspace(MinPfa,MaxPfa,NumPoints);

figure(1)
Pd = privrocpdcalc(Pfa,d,NumPulses,'NonfluctuatingCoherent');

SNRlen = numel(d);
Pd = zeros(numel(Pfa),SNRlen);  % preallocate

erfcinv2PFA = erfcinv(2*Pfa);
for k = 1:SNRlen,
    Pd(:,k) = 0.5*erfc( erfcinv2PFA-(sqrt(d(k)*NumPulses)) );
end
figure(3)
% hold on
semilogx(Pfa,Pd)
grid on
% hold off
%%
% It can be seen from the figure that the measured _Pd_=0.1390 and
% _Pfa_=0.0009 obtained above for the SNR value of 3 dB match a theoretical
% point on the ROC curve.

%% Single Sample Detection Using Noncoherent Receiver
Pfa = 1e-3;


% A noncoherent receiver does not know the phase of the received signal,
% therefore, for target present case, the signal x contains a phase term
% and is defined as

% simulate the signal
x = s.*exp(1i*2*pi*rand(rstream,1,Ntrial)) + n;
y = mf'*x;

%%
% When the noncoherent receiver is used, the quantity used to compare with
% the threshold is the power (or magnitude) of the received signal after
% the matched filter.  In this simulation, we choose the magnitude as the
% sufficient statistic.

z = abs(y);

%%
% Given our choice of the sufficient statistic _z_, the threshold is
% related to _Pfa_ by the equation
%
% $$P_{fa}={\rm exp}\left(-\frac{T^2}{NM}\right)={\rm exp}(-{\rm SNR})\ .$$
%
% The signal to noise ratio threshold SNR for an NP detector can be
% calculated using npwgnthresh as follows:

% snrthreshold = db2pow(npwgnthresh1(Pfa, 1,'noncoherent'));
snrthreshold = (abs(sqrt(gammaincinv(1-Pfa,npulses))))^2;
%%
% The threshold, _T_, is derived from SNR as before
mfgain = mf'*mf;
threshold = sqrt(npower*mfgain*snrthreshold);

%%
% Again, _Pd_ can then be obtained using
Pd = sum(z>threshold)/Ntrial

%%
% Note that this resulting _Pd_ is inferior to the performance we get from
% a coherent receiver.
%
% For the target absent case, the received signal contains only noise.  We
% can calculate the _Pfa_ using Monte-Carlo simulation as

x = n;
y = mf'*x;
z = abs(y);
Pfa = sum(z>threshold)/Ntrial

%%
% The ROC curve for a noncoherent receiver is plotted as
figure(2)
rocsnr1(snrdb,'SignalType','NonfluctuatingNoncoherent','MinPfa',1e-6);

MinPfa = 1e-6;
MaxPfa = 1;
NumPoints =101;
NumPulses=1;

d = db2pow(SNR);
MinPfa = log10(MinPfa);  % convert to log scale for even space sampling
MaxPfa = log10(MaxPfa);
Pfa = logspace(MinPfa,MaxPfa,NumPoints);
SNRlen = numel(SNR);
Pd = zeros(numel(Pfa),SNRlen);  % preallocate
T = real(gammaincinv(1-Pfa,NumPulses));
T = T(:);
SQRT2T = sqrt(2*T);

for k = 1:SNRlen,
    NX = NumPulses*d(k);
    NXT = NX*T;
    Pd(:,k) = marcumq(sqrt(2*NX),SQRT2T);
    Pdtemp = 0;
    if SNR(k)~=0 && SNR(k)~=inf
        for r = 2:NumPulses
            Pdtemp = Pdtemp + (T/NX).^((r-1)/2).*besseli(r-1,2.*sqrt(NXT),1);
        end
    end
    Pd(:,k) = Pd(:,k)+exp(-(T+NX)+2.*sqrt(NXT)).*Pdtemp;
end
figure(4)
% hold on
semilogx(Pfa,Pd)
grid on
% hold off