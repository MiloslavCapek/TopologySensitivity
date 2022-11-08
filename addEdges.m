function [Y, globIndsBF] = addEdges(Y, Z, globIndsBF, globIndsBFAdd)
%% addEdges: add edges specified in globIndsBFAdd back to the structure
%           represented by globIndsBF
% 
% Inputs:
%   Y              ~ actual admittance matrix of a system
%   Z              ~ initial (complete) system (impedance) matrix
%   globIndsBF     ~ actual representation of structure (global positions)
%   globIndsBFAdd  ~ edges to add (global positions)
% 
% Outputs:
%   Y              ~ updated admittance matrix of a system
%   globIndsBF     ~ updated representation of structure (global positions)
% 
% (The code is started from START.m.)
% 
% See also: removeEdges
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

N = length(Z);
% find all edges which are not present in the structure right now
globRemovedEdges = setdiff(1:N, globIndsBF);
% find those edges which have to be added back and their positions
positionsToAdd   = ismember(globRemovedEdges, globIndsBFAdd);
indicesToAdd     = globRemovedEdges(positionsToAdd);

% add edges from the list one by one
for thisEdge = 1:length(indicesToAdd)
    addEdge(indicesToAdd(thisEdge));
end

    function addEdge(m)
        % index-out auxiliary variables
        zn  = Z(globIndsBF, m);
        Znn = Z(m, m);
        
        % block-inversion for fast addition
        k = Y*zn;
        T = 1/(Znn - zn.'*k);
        Y = T*[Y/T + k*k.', -k; ...
                      -k.', 1];
        
        % sort edges in order
        [globIndsBF, retrace] = sort([globIndsBF m]);
        Y = Y(retrace, retrace);
    end
end