function globIndsTR = globBF2globTR(BF, globIndsBF)
%% globBF2globTR: find triangles (global indices) adjacent to basis
%                functions (global indices)
% 
% Inputs:
%   BF             ~ basis function structure (from AToM, see [1])
%   globIndsBF     ~ actual representation of structure (global positions)
% 
% Outputs:
%   globIndsTR     ~ actual representation of structure in terms of
%                    triangles
% 
% (The code is started from START.m.)
% 
% [1] AToM: Antenna Toolbox for MATLAB, [on-line]: www.antennatoolbox.com
% 
% See also: globTR2globBF
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

globIndsTR = unique(BF.data(globIndsBF, [2 4])); % [trPlus trMinus]