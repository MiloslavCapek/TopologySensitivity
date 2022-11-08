function thisLocIndsBF = glob2loc(globIndsBF, thisGlobIndsBF)
%% glob2loc: recalculates global positions (indices) to local (those
%           actually used)
% 
% Inputs:
%   globIndsBF     ~ actual representation of structure (global positions)
%   thisGlobIndsBF ~ a vector of global indices
% 
% Outputs:
%   thisLocIndsBF ~ corresponding vector of local indices
% 
% (The code is started from START.m.)
% 
% See also: loc2glob
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

[~, ~, thisLocIndsBF] = intersect(thisGlobIndsBF, globIndsBF, 'stable');
thisLocIndsBF = thisLocIndsBF.';