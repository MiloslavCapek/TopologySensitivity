function P = evaluateQuadraticForm(A, I, globIndsBF, globIndsBFExtra)
%% evaluateLinearForm: evaluate linear form (A*I) for all currents
% 
% Inputs:
%   A               ~ matrix operator to be used for I'*A*I
%   I               ~ current vectors (one or many, same but arbitrary
%                     length)
%   globIndsBF      ~ actual representation of structure (global positions)
%   globIndsBFExtra ~ indexes the I entries into A matrix
%                     (optional for removal, mandatory for addition)
% 
% Outputs:
%   P               ~ column of values
% 
% (The code is started from START.m.)
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

% use different strategies (for different cost) based on complexity...
nInputs = nargin;
if nInputs < 3 % for topo. sens. with edge removal, A, I have the same size
    P = sum(dot(I, A*I), 1).';
elseif nInputs < 4 % for topo. sens. with edge removal, A has original size
    P = sum(dot(I, A(globIndsBF, globIndsBF)*I), 1).';
else % for topology senstivity with edge addition (can be accelerated?!?)
    B  = size(I, 2);
    N  = size(A, 2);
    I_ = zeros(N, B);
    for b = 1:B
        I_([globIndsBF.'; globIndsBFExtra(b)], b) = I(:, b);
    end
    
% reduce size of matrices
    inds = setdiff(1:N, [globIndsBF, reshape(globIndsBFExtra, 1, [])]);
    I_(inds, :) = [];
    A(inds, :) = [];
    A(:, inds) = [];
    
% evaluate (potentially reduced) quadratic form
    P = sum(dot(I_, A*I_), 1).';
end