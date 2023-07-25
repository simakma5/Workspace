function S = transmLineSPar(gamma, L, Zc, Z0)
%% transmLineSPar computes S-parameters of segment of transmission line with
% propagation constant gamma, length L, impedance Zc measured with ports with
% impedance Z0.
% NOTE: third and upper dimension of resulting variable S can have whatever
% lengths, a computations is fully vectorized
% Implemented according to User’s Guide to S-Parameters Explorer 2.0
% (http://www.eecircle.com/downloads/SparViewer/SPEX20_guide.pdf) p. 31
%
%  INPUTS
%   gamma: propagation constant of transmission line in form alpha + j*beta,
%          alpha-Attenuation const., beta-Phase const., i.e. S21 of Line 
%          (measured with ports with impedance Zc) is exp(-gamma*L), 
%          double [1 x 1 x nFreq]
%   L: length of transmission line in meters, double [1 x 1]
%   Zc: impedance of transmission line, double [1 x 1] in case of
%       non-disspersive line, double [1 x 1 x nFreq] in dispersive case
%   Z0: impedance of ports which are utilized for measuremnts of transmission
%       line, double [1 x 1]
%
%  OUTPUTS
%   S: S-parameters of segment of transmission line, double [2 x 2 x nFreq]
%
% © 2021, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz

theta = Zc./Z0;

gammaL = gamma.*L;
S11 = (sinh(gammaL).*(theta - 1./theta)) ...
   ./(2*cosh(gammaL) + sinh(gammaL).*(theta + 1./theta));
S21 = 2./(2*cosh(gammaL) + sinh(gammaL).*(theta + 1./theta));

S = [S11, S21; S21, S11];

end