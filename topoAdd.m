function [IB, globIndsBF_IB] = topoAdd(...
    Z, V, Y, I, globIndsBF, globIndsBFTopo)
%% topoAdd: evaluate all possible topology changes subject to addition
% 
% Inputs:
%   Z              ~ initial (complete) system (impedance) matrix
%   V              ~ initial (complete) excitation vector
%   Y              ~ actual admittance matrix of a system
%   I              ~ actual current, i.e., I = Y*V(globIndsBF)
%   globIndsBF     ~ actual representation of structure (global positions)
%   globIndsBFTopo ~ edges to investigate (global positions)
% 
% Outputs:
%   IB             ~ matrix of column vectors containing currents, one by
%                    one corresponding to particular topology change
%   globIndsBF_IB  ~ column vectors of global indices corresponding to each
%                    topology change (column by column)
% 
% (The code is started from START.m.)
% 
% See also: topoRemove
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

% if the last input is not specified, evaluate all possibilities
if nargin < 6
    globIndsBFTopo = setdiff(1:size(Z, 1), globIndsBF);
end

% linear index into diagonal
ZdiagIndsBF = size(Z,1)*(globIndsBFTopo-1) + globIndsBFTopo;

% index-out auxiliary variables
ZB  = Z(globIndsBF, globIndsBFTopo);
X   = Y*ZB;

% complicated version optimized for comp. time (formulas are in the paper)
vm = (X.'*V(globIndsBF) - V(globIndsBFTopo)).';
IB = [I; 0] + (vm ./ (Z(ZdiagIndsBF) - sum(ZB .* X, 1))) .* ... vector1
    [X; -1*ones(size(globIndsBFTopo))];

% if the indices are desired by the user as well, compute them
if nargout > 1
    globIndsBF_IB = ...
        [repmat(globIndsBF.', [1, length(globIndsBFTopo)]); ...
         globIndsBFTopo];
end