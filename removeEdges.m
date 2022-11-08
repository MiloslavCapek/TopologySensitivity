function [Y, globIndsBF] = removeEdges(Y, globIndsBF, globIndsBFRem)
% removeEdges: remove edges specified in globIndsBFRem back to the 
%              structure represented by globIndsBF
% 
% Inputs:
%   Y              ~ actual admittance matrix of a system
%   globIndsBF     ~ actual representation of structure (global positions)
%   globIndsBFRem  ~ edges to remove (global positions)
% 
% Outputs:
%   Y              ~ updated admittance matrix of a system
%   globIndsBF     ~ updated representation of structure (global positions)
% 
% (The code is started from START.m.)
% 
% See also: addEdges
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

% if the last input is not specified, evaluate all possibilities
indicesBFToRemove = find(ismember(globIndsBF, globIndsBFRem));

% remove edges from the list one by one
for thisEdge = 1:length(indicesBFToRemove)
    removeEdge(indicesBFToRemove(thisEdge));
end

% update admittance matrix
Y(indicesBFToRemove, :) = [];
Y(:, indicesBFToRemove) = [];

% update representation of the structure
globIndsBF(indicesBFToRemove) = [];

    % application of Sherman-Morrison-Woodbury (can be vectorized)
    function removeEdge(m)
        % division by Y(m,m) is advantageously done for vector
        Y = Y - (Y(:,m)/Y(m,m)) * Y(:,m).';
    end
end