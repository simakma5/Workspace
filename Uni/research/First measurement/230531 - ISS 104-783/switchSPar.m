function SSwitch = switchSPar(S)
%% switchSPar change direction of connection of S-parameter box.
% In 2-Port case it switches ports 1-2.
% In 3-Port case it change ports 1-3, 2-2.
% In 4-Port case it schitches ports 1-4, 2-3.
%
%  INPUTS
%   S: S-parameters of circuit, double [nPorts x nPorts x nFreq]
%
%  OUTPUTS
%   SSwitch: switched S-parameters, double [nPorts x nPorts x nFreq]
%
% © 2020, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz

SSwitch = rot90(S, 2);

end