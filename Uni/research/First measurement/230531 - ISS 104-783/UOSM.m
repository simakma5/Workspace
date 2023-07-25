function [Ea, Eb, maxPhaseThruDev] = UOSM(P1Meas, P2Meas, SThruMeas, P1Model, P2Model, ...
   thruS21Model, freq)
%% UOSM compute error model of 2-port VNA.
% Error boxes A and B are computed from OSM method. Additional transmission
% normalization of error box B is made in accordance to Thru measurement. It is
% assumed that port n. 1 of box A is connected to VNA, port n. 2 is connected to
% DUT. Port n. 1 of box B is connected to DUT, port n. 2 is connected to VNA.
% This use 7-term error model. Error box A can be reciprocal or not depending on
% OSM calibration implementation. Reflection coeffitients of error boxes
% computed from OSM calibration are not changed. Transmission coeficients of box
% B are somehow changed.
% Assumed phase of S21 of Thru is utilized to distinguish correct solution of
% alpha.
% Thru standard is assumed to be reciprocal (S21=S12). It can have some
% reflections, but ussually small because of linearity of phase of transmission.
% It is necessary to ensure frequency step providing smaller phase step than
% pi/2. Minimal number of frequency points is two.
% This model is not changed even when the load impedances inside VNA are changed
% because of switching from reveiver to source, i.e. on contrary to 12-term
% model it has no separate definition for forward and reverse measurement.
%
%  INPUTS
%   P1Meas, P2Meas: struct with fields:
%         .open: measurement of Open standard on port 1 and 2, double [nFreq x 1]
%         .short: measurement of Short standard on port 1 and 2, double [nFreq x 1]
%         .match: measurement of Match standard on port 1 and 2, double [nFreq x 1]
%   thruMeas: measurement of Thru standard, double [2 x 2 x nFreq]
%   P1Model, P2Model: struct with fields:
%          .open: model of Open standard on port 1 and 2, double [nFreq x 1]
%          .short: model of Short standard on port 1 and 2, double [nFreq x 1]
%   thruS21Model: estimation of Thru standard transmission, it is necessary to ensure
%                 frequency step providing smaller phase step than pi/2, 
%                 double [nFreq x 1]
%
%  OUTPUTS
%   Ea, Eb: error boxes, double [2 x 2 x nFreq]
%   maxPhaseThruDev: deviation of phase between thruS21Model and calibrated Thru
%                    in degrees, double [nFreq x 1]
%
% Implemented according to: Ferrero, A., Pisani, U.: Two-Port Network Analyzer
%                           Calibration Using an Unknown "Thru".
%
% © 2017, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz

openP1Meas = P1Meas.open;
shortP1Meas = P1Meas.short;
matchP1Meas = P1Meas.match;

openP2Meas = P2Meas.open;
shortP2Meas = P2Meas.short;
matchP2Meas = P2Meas.match;

openP1Model = P1Model.open;
shortP1Model = P1Model.short;
matchP1Model = P1Model.match;

openP2Model = P2Model.open;
shortP2Model = P2Model.short;
matchP2Model = P2Model.match;

nFreq = length(openP1Meas);
%% OSM calibration

Ea = RRR1PortCalibration(openP1Meas, shortP1Meas, matchP1Meas, ...
   openP1Model, shortP1Model, matchP1Model, freq, 'n');
Eb = RRR1PortCalibration(openP2Meas, shortP2Meas, matchP2Meas, ...
   openP2Model, shortP2Model, matchP2Model, freq, 'n');
%% Get correct alpha solution from UOSM method

Ya = s2t(Ea)./Ea(1, 2, :); % Ok
Yb = switchSPar(s2t(Eb)./Eb(1, 2, :)); % Ok

TThruMeas = s2t(SThruMeas);
X = zeros(size(SThruMeas), 'like', 1j);
detTThruMeas = zeros(nFreq, 1, 'like', 1j);
detYa = zeros(nFreq, 1, 'like', 1j);
detYb = zeros(nFreq, 1, 'like', 1j);

% compute determinants and theoretical Thru parameters utilizing just OSM
% calibration results
for iFreq = 1:nFreq
   detTThruMeas(iFreq) = det(TThruMeas(:, :, iFreq));
   detYa(iFreq) = det(Ya(:, :, iFreq));
   detYb(iFreq) = det(Yb(:, :, iFreq));
   X(:, :, iFreq) = Ya(:, :, iFreq)\TThruMeas(:, :, iFreq)*Yb(:, :, iFreq);
end
alpha1 = sqrt(detTThruMeas.*detYb./detYa);
alpha2 = -alpha1;
% get two solutions of S21 of Thru standard
S21ThruEst1 = alpha1./squeeze(X(2, 2, :));
S21ThruEst2 = alpha2./squeeze(X(2, 2, :));

logIndSol1 = abs(angle(S21ThruEst1./thruS21Model)) <= abs(angle(S21ThruEst2./thruS21Model));
alpha = alpha2;
alpha(logIndSol1) = alpha1(logIndSol1);
maxPhaseThruDev = (angle(thruS21Model./(alpha./squeeze(X(2, 2, :))))/pi*180);

for iFreq = 1:nFreq
   Yb(:, :, iFreq) = alpha(iFreq)*inv(Yb(:, :, iFreq)); % inv is ok
end

Eb = t2s(Yb);


% % get steps in phase of S21 of Thru (should be smooth)
% dS21ThruPhase1 = angle(S21ThruEst1(2:end)./S21ThruEst1(1:end-1));
% % dS21ThruPhase2 = angle(S21ThruEst2(2:end)./S21ThruEst2(1:end-1));
% 
% % find huge steps in phase larger than pi/2
% % on places where alpha1 substitutes alpha2 is step pi
% logIndHugePhaseStep = abs(dS21ThruPhase1) >= pi/2;
% indUniqueSol1 = logical(mod(cumsum(logIndHugePhaseStep), 2));
% indUniqueSol1 = [indUniqueSol1(1); indUniqueSol1];
% indUniqueSol2 = ~indUniqueSol1;
% 
% % sort solution 1 and 2
% alphaSol1 = zeros(size(alpha1));
% alphaSol1(indUniqueSol1) = alpha1(indUniqueSol1);
% alphaSol1(indUniqueSol2) = alpha2(indUniqueSol2);
% 
% alphaSol2 = zeros(size(alpha1));
% alphaSol2(indUniqueSol2) = alpha1(indUniqueSol2);
% alphaSol2(indUniqueSol1) = alpha2(indUniqueSol1);
% 
% % two solutions for Eb
% Tb1 = shiftdim(alphaSol1, -2).*Yb;
% Tb2 = shiftdim(alphaSol2, -2).*Yb;
% 
% % get two estimations of S21 of thru
% S21PhaseThruSol1 = alphaSol1./squeeze(X(2, 2, :));
% S21PhaseThruSol2 = alphaSol2./squeeze(X(2, 2, :));
% 
% % check mean deviation of phase between estimated and expected S21
% dS21PhaseSol1 = mean(abs(angle(S21PhaseThruSol1./thruS21Model)));
% dS21PhaseSol2 = mean(abs(angle(S21PhaseThruSol2./thruS21Model)));
% 
% % get solution with smaller deviation
% if dS21PhaseSol1 < dS21PhaseSol2
%    Eb = t2s(Tb1);
% else
%    Eb = t2s(Tb2);
% end

end