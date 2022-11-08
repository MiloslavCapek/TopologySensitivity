function [globIndsBFRem, globIndsBFAdd] = objectBoundary(BF, globIndsBF)
%% objectBoundary: find boundary of an object and corresponding basis 
%                  functions which might be subject of addition/removal
%                  along the boundary
% 
% Inputs:
%   BF             ~ basis function structure (from AToM, see [1])
%   globIndsBF     ~ actual representation of structure (global positions)
% 
% Outputs:
%   globIndsBFRem  ~ global indices of edges (basis functions) which might
%                    be removed to get the topology senstivity along the
%                    boundary of an object
%   globIndsBFAdd  ~ the same as "globIndsBFRem" but for additions
% 
% This function can be applied to significantly accelerate the optimization
% by optimizing only along the boundary of an object.
% 
% (The code is started from START.m.)
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

globIndsBFout = setdiff(1:BF.nUnknowns, globIndsBF);

globIndsTRin  = globBF2globTR(BF, globIndsBF).';
globIndsTRout = globBF2globTR(BF, globIndsBFout);

globIndsTRboundary = intersect(globIndsTRin, globIndsTRout);

globIndsBFall = ...
    BF.data(any(ismember(BF.data(:, [2 4]), globIndsTRboundary), 2), 7).';

globIndsBFRem = intersect(globIndsBFall, globIndsBF);
globIndsBFAdd = intersect(globIndsBFall, globIndsBFout);