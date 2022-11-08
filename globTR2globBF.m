function globIndsBF = globTR2globBF(BF, globIndsTR)
%% globTR2globBF: find triangles (global indices) adjacent to basis
%                functions (global indices)
% 
% Inputs:
%   BF             ~ basis function structure (from AToM, see [1])
%   globIndsTR     ~ actual representation of structure in terms of
%                    triangles
% 
% Outputs:
%   globIndsBF     ~ actual representation of structure (global positions)
% 
% (The code is started from START.m.)
% 
% [1] AToM: Antenna Toolbox for MATLAB, [on-line]: www.antennatoolbox.com
% 
% See also: globBF2globTR
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

globIndsBF = BF.data(all(ismember(BF.data(:, [2 4]), globIndsTR), 2), 7).';