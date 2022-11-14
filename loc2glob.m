function thisGlobIndsBF = loc2glob(globIndsBF, thisLocIndsBF)
%% loc2glob: recalculates local positions (indices) to global (those actually
%           used)
% 
% Inputs:
%   globIndsBF     ~ actual representation of structure (global positions)
%   thisLocIndsBF ~ corresponding vector of local indices
% 
% Outputs:
%   thisGlobIndsBF ~ a vector of global indices
% 
% (The code is started from START.m.)
% 
% See also: glob2loc
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

thisGlobIndsBF = globIndsBF(thisLocIndsBF);