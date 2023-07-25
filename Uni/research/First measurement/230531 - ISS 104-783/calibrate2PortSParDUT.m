function S_DUT = calibrate2PortSParDUT(Sa, Sb, S_DUT_m)
%% calibrate2PortSParDUT use two error boxes Sa and Sb to compute S-parameters
% S_DUT from measured S_DUT_m.
% The 8-term error model is used. It can apply n error boxes to one DUT. Port n. 1
% of error box A is connected to port 1 of VNA, port n. 2 to port n. 1 of DUT. Port n. 1
% of error box B is connected to port 2 of DUT and port 2 of B is connected to
% port 2 of VNA. It is necessary the DUT to have non-zero transmission S21.
%
%  INPUTS
%   Sa, Sb: S-parameters of error boxes, double [2 x 2 x nFreq x n]
%   S_DUT_m: measured S-parameters (Sa - DUT - Sb),
%            double [2 x 2 x nFreq]
%
%  OUTPUTS
%   S_DUT: S-parameters of calibrated 2-port DUT, double [2 x 2 x nFreq x n]
%
% © 2017, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz


%% Calibrating DUT

nCal = size(Sa, 4);
S_DUT = zeros(size(Sa));
T_DUT_m = s2t(S_DUT_m);
nFreq = size(S_DUT_m, 3);

for iCal = 1:nCal
   Ta = s2t(Sa(:, :, :, iCal));
   Tb = s2t(Sb(:, :, :, iCal));
   T_DUT = zeros(size(S_DUT_m));
   
   for iFreq = 1:nFreq
      %    T_DUT_cal(:, :, iFreq) = inv(Ta(:, :, iFreq))*T_DUT_m(:, :, iFreq)*inv(Tb(:, :, iFreq));
      T_DUT(:, :, iFreq, iCal) = Ta(:, :, iFreq)\T_DUT_m(:, :, iFreq)/Tb(:, :, iFreq);
   end
   
   S_DUT(:, :, :, iCal) = t2s(T_DUT(:, :, :, iCal));
   
end

S_DUT = squeeze(S_DUT);

end

function t = s2t(s)
t = zeros(size(s));
t(1,1,:) = s(1, 2,:) - s(1, 1,:).*s(2,2,:)./s(2,1,:);
t(1,2,:) = s(1, 1,:)./s(2,1,:);
t(2,1,:) = -s(2,2,:)./s(2,1,:);
t(2,2,:) = 1./s(2,1,:);
end

function s = t2s(t)
s = zeros(size(t));
s(1,1,:) = t(1,2,:)./t(2,2,:);
s(1,2,:) = t(1,1,:) - t(1,2,:).*t(2,1,:) ./t(2,2,:);
s(2,1,:) = 1./t(2,2,:);
s(2,2,:) = -t(2,1,:) ./t(2,2,:);
end