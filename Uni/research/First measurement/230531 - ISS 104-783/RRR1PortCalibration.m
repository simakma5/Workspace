function Err = RRR1PortCalibration(R1Meas, R2Meas, R3Meas, R1Model, R2Model, ...
   R3Model, freq, recOrNonRec)
%% RRR1PortCalibration computes error two-port Err of 1-port VNA.
% As calibratios standards are used 3 Reflects. It is not specified what kind of
% standard can be used. On port n. 1 of Err VNA is connected, on port n. 2 DUT
% (caliber) is connected. Implemented according to Handbook of Microwave
% Component Measurements, Dunsmore, 2012, Wiley, p. 132.
%
%  INPUTS
%   R1Meas, R2Meas, R3Meas: measured S-parameters of 3 Reflect standards, 
%                           double [nFreq x 1]
%   R1Model, R2Model, R3Model: true S-parameters of Reflect standards in 
%                              reference plane, double [nFreq x 1]
%   recOrNonRec: 'r' for reciprocal circuit, i.e. E12 and E21 are the same, 
%                'n' for non-reciprocal circuit, where E21 = E21*E12, char [1 x 1]
%                optional parameter, default is 'r'
%
%  OUTPUTS
%   Err: S-parameters of error two-port, double [2 x 2 x nFreq]
%
% © 2019, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz


%% computation

nFreq = length(freq);
Err = zeros(2, 2, nFreq, 'like', 1j);
ERF = zeros(size(freq), 'like', 1j);
for iFreq = 1:nFreq
   GammaAR1 = R1Model(iFreq);
   GammaAR2 = R2Model(iFreq);
   GammaAR3 = R3Model(iFreq);
   
   GammaMR1 = R1Meas(iFreq);
   GammaMR2 = R2Meas(iFreq);
   GammaMR3 = R3Meas(iFreq);
   
   A = [1 GammaAR1 GammaAR1*GammaMR1;
        1 GammaAR2 GammaAR2*GammaMR2;
        1 GammaAR3 GammaAR3*GammaMR3];
   B = [GammaMR1; GammaMR2; GammaMR3];
   x = A\B;
   
   Err(1, 1, iFreq) = x(1); % EDF
   Err(2, 2, iFreq) = x(3); % ESF
   ERF(iFreq) = x(2) + x(1)*x(3);
end

if ~exist('recOrNonRec', 'var')
   recOrNonRec = 'r';
end

switch recOrNonRec
   case 'r'
      % make S21 and S12 of error box the same
      Err(2, 1, :) = sqrtOfComplexTransmission(ERF, freq);
      Err(1, 2, :) = Err(2, 1, :);
   case 'n'
      Err(2, 1, :) = ERF;
      Err(1, 2, :) = ones(size(freq));
   otherwise
      error('Wrong reciprocity input parameter.');
end

end