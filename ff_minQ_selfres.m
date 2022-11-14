function FF = ff_minQ_selfres(OP, I, globIndsBF, globIndsBF_tested, BFtype)
%% ff_minQ_selfres: evaluates antenna Q factor with penalty for
%                  non-self-resonance currents
% 
% This function evaluates antenna Q factor and adds penalty term
% 
%      FF = ((max(omWM, omWE) + alpha*abs(omWM - omWE)) ./ Pr),
% 
% where omWM is stored magnetic, omWE stored electric energies (both
% multiplied by angular frequency), Pr is radiated power, and alpha is
% user-defined weight (higher value means higher preference of self-resonance).
% 
% Inputs:
%   OP                ~ MATLAB structure containing all required 
%                       information about the structure. Here, matrices 
%                       omXm, omXe, Z, and S1 are required. 
%                       S1 is a factorization of radiation part of Z, i.e., 
%                       S1'*S1 = real(Z) for lossless structures.
%   I                 ~ actual currents, matrix of arbitrarily many columns
%   globIndsBF        ~ actual representation of structure (global positions)
%   globIndsBF_tested ~ edges which are subject of evaluation (empty / a
%                       vector)
%   BFtype            ~ type of currents used for evaluation of the fitness
%                       function, 0 = no topology modifications, -1  = test 
%                       on removals, +1 = test on additions
% 
% Outputs:
%   FF                ~ fitness function values (in a column vector)
% 
% (The code is started from START.m.)
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

%% (Quadratic/Linear-based quantities)
% Code might seem complicated, but implemented for high versatility (same
% functions for arbitrary metric) and high speed:
nBFgiven = size(globIndsBF, 2);
if BFtype == -1 || BFtype == 0 % Investigation of edge removal
    omWM = real(evaluateQuadraticForm(OP.omXm, I, globIndsBF));
    omWE = real(evaluateQuadraticForm(OP.omXe, I, globIndsBF));
    if nBFgiven <= size(OP.S,1) % I'*A*I
        Pr = real(evaluateQuadraticForm(real(OP.Z), I, globIndsBF));
    else % |B*I|.^2, where B*B' = A (low-rank approximation)
        SI = evaluateLinearForm(OP.S, I, globIndsBF);
        Pr = dot(SI, SI, 1).';
    end
elseif BFtype == +1 % Investigation of edge addition
    nBFtested = size(globIndsBF_tested, 2);
    omWM = real(evaluateQuadraticForm(OP.omXm, I, globIndsBF, globIndsBF_tested));
    omWE = real(evaluateQuadraticForm(OP.omXe, I, globIndsBF, globIndsBF_tested));
    if nBFgiven+nBFtested <= size(OP.S,1) % I'*A*I
        Pr = real(evaluateQuadraticForm(real(OP.Z), I, globIndsBF, globIndsBF_tested));
    else % |B*I|.^2, where B*B' = A (low-rank approximation)
        SI = evaluateLinearForm(OP.S, I, globIndsBF, globIndsBF_tested);
        Pr = dot(SI, SI, 1).'; % Pr  = sum(abs(SI).^2, 1).';
    end
end

% fitness function
alpha = 1.5; % penalization to reach self-resonance
FF = ((max(omWM, omWE) + alpha*abs(omWM - omWE)) ./ Pr);

end