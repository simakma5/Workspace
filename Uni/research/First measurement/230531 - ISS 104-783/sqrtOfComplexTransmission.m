function [sqrS21, phaseErrDeg] = sqrtOfComplexTransmission(S21, freq, maxExtrapFreq)
%% sqrtOfComplexTransmission compute square root of complex transmission.
% Square root computation of complex transmission is usual when implementing VNA
% calibration methods. It is assumed that the phase of transmission coefficient
% is evolving linearly (non-dispersive case). To obtain a unambibuous solution
% it is necessary extrapolate unwrapped phase to DC and then make correction of
% the phase to represent a specific delay circuit.
% NOTE: works only for non-disspersive transmissions!
% NOTE: transmission has to have phase step less than 180 deg. per freq. point.
% TODO: enable extrapolation of disspersive transmissions (waveguides).
% TODO: compute the estimated delay and return it as second parameter.
%
%  INPUTS
%   S21: transmission coefficient, complex double [nFreq x 1]
%   freq: frequency points of transmission, double [nFreq x 1]
%   maxExtrapFreq: optional input, maximal frequency of S21 samples which are 
%                  used for extrapolation, default: freq(end), double [1 x 1]
%
%  OUTPUTS
%   sqrS21: sqrt of transmission, double [nFreq x 1]
%   phaseErrDeg: phase at zero frequency observed by linear extrapolation,
%                theoretically should be zero, double [1 x 1]
%
% © 2019, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz


%% inputs check
if ~exist('maxExtrapFreq', 'var')
   % utilize all samples of S21
   nFreq = length(freq);
   indForExtrap = 1:nFreq;
else
   [~, I] = min(abs(freq - maxExtrapFreq));
   indForExtrap = 1:I;
end

%% computation
S21Phase = angle(S21);
unwrappedPhase = unwrap(S21Phase);
% extrapolate phase to DC
% it is actually the phase offset caused by wrapping the phase from the non-zero
% frequency
orderOfApprox = 1; % 1 for linear extrapolation
V = freq(indForExtrap).^(orderOfApprox:-1:0);
p = V\unwrappedPhase(indForExtrap);
phaseOffset = polyval(p, 0); % extrapolated phase at DC (zero frequency)
% check for the interpolated value: vals = polyval(p, freq);
% solution with interp1 use just first two points to compute the extrapolation
% to zero (no least square solution)
% phaseOffset = interp1(freq, unwrappedPhase, 0, 'linear', 'extrap');
% round offset to nearest pi
n2Pi = round(phaseOffset/(2*pi));
% the phase should start exactly with zero phase
% in practice there is always some deviation
% correction of the phase
phaseCorr = n2Pi*2*pi;
phaseErrDeg = (phaseOffset - phaseCorr)/pi*180;
% show warning for excesively large deviation
% if phaseErrDeg >= 10
%    warning('Extrapolation of phase of transmission to DC has deviation of %.1f deg.\n', ...
%       phaseErrDeg);
% end
% make correction
corrUnwrPhase = unwrappedPhase - phaseCorr;
% make square root computation as division of the unwrapped phase
sqrS21 = sqrt(abs(S21)).*exp(1j*corrUnwrPhase./2);

end