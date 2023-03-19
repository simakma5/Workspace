function w = CPWW(substrate, GND, g, f0, Z, dZ, w1, w2, dw)
%% CPWW computes width of a strip of coplanar waveguide by optimization of txlineCPW.
% It requres RF Toolbox.
%
%  INPUTS
%   substrate: definition of the substrate parameters as a cell array in form 
%              of comma-separated pairs {'name1', value1, 'name2', value2}, name
%              of parameters are EpsilonR (-), Height (m), LossTangent (-), 
%              Thickness (m) and SigmaCond (S/m), cell [1 x 2*nPar]
%   GND: presence of infinite-bottom ground plane, logical [1 x 1]
%   g: width of the slot between strip and groupd plates in meters, double [1 x 1]
%   f0: frequency in Hz where the impedance holds, double [1 x 1]
%   Z: wanted characteristic impedance of the CPW line, double [1 x 1]
%   dZ: precission of the reached impedance as stop criteria of the optimizer, 
%       e.g. 0.1 Ohm, double [1 x 1]
%   w1, w2: bounds of the width of the microstrip line, w1 < w2, double [1 x 1]
%   dw: minimal step in width of the microstrip line as stop criteria of the 
%       optimizer, e.g. 1e-6, double [1 x 1]
%
%  OUTPUTS
%   w: width of the strip in coplanar waveguide with impedance Z, double [1 x 1]
%
%  EXAMPLE
%   % find width of a strip in CPWG line with impedance 50 Ohm fabricated on
%   % RO4350B substrate with precision of 0.1 Ohm in range from 0.1 to 7 mm
%   substrate = {'EpsilonR', 3.66, 'Height', 0.508e-3, 'LossTangent', 0.0037, ...
%                'Thickness', 17e-6, 'SigmaCond', 5.88e7};
%   w = CPWW(substrate, true, 0.2e-3, 10e9, 50, 0.1, 0.1e-3, 7e-3, 1e-6)
%
% © 2022, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz

w = fminbnd(@(w)fitnessFcn(substrate, GND, g, f0, Z, w), w1, w2, ...
   optimset('TolFun', dZ, 'TolX', dw));

end

function f = fitnessFcn(substrate, GND, g, f0, Z, w)
line = txlineCPW(substrate{:}, 'ConductorBacked', GND, 'SlotWidth', g, ...
   'ConductorWidth', w);
ZLine = getZ0(line, f0);
f = abs(ZLine - Z);
end