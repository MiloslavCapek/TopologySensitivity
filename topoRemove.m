function IB = topoRemove(Y, I, globIndsBF, globIndsBFTopo)
%% topoRemove: evaluate all possible topology changes subject to removals
% 
% Inputs:
%   Y              ~ actual admittance matrix of a system
%   I              ~ actual current, i.e., I = Y*V(globIndsBF)
%   globIndsBF     ~ actual representation of structure (global positions)
%   globIndsBFTopo ~ edges to investigate (global positions)
% 
% Outputs:
%   IB             ~ matrix of column vectors containing currents, one by
%                    one corresponding to particular topology change
% 
% (The code is started from START.m.)
% 
% See also: topoAdd
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

% no transpose here (.') since column vector realized by multiplication
Ydiag = diag(Y);

% use following lines depending on how many changes are investigated
if nargin < 4 % for all (relatively faster code)
    IB = I - Y .* (I ./ Ydiag).';
else % for a particular list
    locIndsTopo = glob2loc(globIndsBF, globIndsBFTopo);
    IB = I - Y(:,locIndsTopo) .* (I(locIndsTopo)./Ydiag(locIndsTopo)).';
end