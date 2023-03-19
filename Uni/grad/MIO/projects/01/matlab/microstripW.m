function w = microstripW(substrate, f0, Z, dZ, w1, w2, dw)
%% microstripW computes width of a microstrip line by optimization of txlineMicrostrip.
% It requres RF Toolbox.
%
%  INPUTS
%   substrate: definition of the substrate parameters as a cell array in form 
%              of comma-separated pairs {'name1', value1, 'name2', value2}, name
%              of parameters are EpsilonR (-), Height (m), LossTangent (-), 
%              Thickness (m) and SigmaCond (S/m), cell [1 x 2*nPar]
%   f0: frequency in Hz where the impedance holds, double [1 x 1]
%   Z: wanted characteristic impedance of the microstrip line, double [1 x 1]
%   dZ: precission of the reached impedance as stop criteria of the optimizer, 
%       e.g. 0.1 Ohm, optional parameter, default 0.01 Ohm, double [1 x 1]
%   w1, w2: bounds of the width of the microstrip line, w1 < w2, optional 
% %         parameter, default w1=0.05*H, w2 = 20*H, double [1 x 1]
%   dw: minimal step in width of the microstrip line as stop criteria of the 
%       optimizer, e.g. 1e-6, optional parameter, default dw = (w2 - w1)/1e3, 
%       double [1 x 1]
%
%  OUTPUTS
%   w: width of microstrip line with impedance Z, double [1 x 1]
%
%  EXAMPLE
%   % find width of a microstrip line with impedance 50 Ohm fabricated on RO4350B 
%   % substrate with precision of 0.1 Ohm in range from 0.1 to 7 mm
%   substrate = {'EpsilonR', 3.66, 'Height', 0.508e-3, 'LossTangent', 0.0037, ...
%                'Thickness', 17e-6, 'SigmaCond', 5.88e7};
%   w = microstripW(substrate, 10e9, 50, 0.1, 0.1e-3, 7e-3, 1e-6)
%   % without defining Z and w limits
%   w = microstripW(substrate, 10e9, 50)
%
% © 2022, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz

if nargin == 3
   dZ = 0.01;
end

if nargin <= 4 % w1, w2, dw is not defined
   % find substrate height index
   substrHeightInd = find(cellfun(@(subProp)strcmp(subProp, 'Height'), substrate)) + 1;
   substrateH = substrate{substrHeightInd};
   % search solution within txlineMicrostrip validity region (w/H in <0.05; 20>)
   w1 = 0.05*substrateH;
   w2 = 20*substrateH;
   dw = (w2 - w1)/1e3;
end

w = fminbnd(@(w)fitnessFcn(substrate, f0, Z, w), w1, w2, ...
   optimset('TolFun', dZ, 'TolX', dw));

end

function f = fitnessFcn(substrate, f0, Z, w)
line = txlineMicrostrip(substrate{:}, 'Width', w);
ZLine = getZ0(line, f0);
f = abs(ZLine - Z);
end