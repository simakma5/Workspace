function SParInCell = convert4DSparToCell(SPar)
%% convert4DSparToCell converts S-parameters from several measurements (in 4th
% dimension) to cell variable usually for usage in easy ploting.
%
%  INOUTS
%   SPar: S-parameters, double [nP x nP x nFreq x nPar]
%
%  OUTPUTS
%   SParInCell: S-parameters arranged in cell, cell [nPar x 1], in every cell is
%               double [nP x nP x nFreq]
%
% © 2020, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz

nP = size(SPar, 1); % the same as dim. 2
nFreq = size(SPar, 3);
nPar = size(SPar, 4);

SParInCell = mat2cell(SPar, nP, nP, nFreq, ones(nPar, 1));
SParInCell = SParInCell(:);

end