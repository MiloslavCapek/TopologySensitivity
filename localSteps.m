function [FF, GenusInit, GenusFinal, GenusHistory, tt, candidates] = ...
    localSteps(OP, optData)
%% localSteps: perform all local steps towards local minimum starting from 
%              a random seed (if called multiple time, it results in a Monte 
%              Carlo analysis). Each local step is performed in greedy sense.
% 
% Inputs:
%   OP           ~ MATLAB structure containing all variables and fields
%                  fully describing the optimization region and all
%                  necessary MoM matrices. Majority of the fields are
%                  evaluated by AToM, see [1]. Mandatory fields are:
%                  OP.Mesh (discretization), OP.BF (basis functions),
%                  OP.ports (position of discrete ports, if any), OP.V
%                  (excitation vector), OP.Z (impedance matrix), 
%                  OP.Zsystem (system matrix - might be equal to OP.Z),
%                  etc. See START.m for details.
%   optData      ~ MATLAB structure containing optimization settings. 
%                  See START.m for details.
% 
% Outputs:
%   FF           ~ row vector of fitness function value in each iteration
%   GenusInit    ~ a gene representing the initial structure
%   GenusFinal   ~ a gene representing the resulting structure
%   GenusHistory ~ a row vector containing history of all structural
%                  changes. positive value: an edge of that global position
%                  was added, negative value: an edge of that global
%                  position was removed
%   tt           ~ a row vector of computational time (commulative)
%   candidates   ~ a vector 1x2 containing numbers of antenna samples 
%                  evaluated, (1,1): investigated for removal, (1,2):
%                  investigated for addition
% 
% (The code is started from START.m.)
% 
% See also: START
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

% Optimization settings:
fitnessFunction = optData.fitness;
protectedEdges  = optData.protectedEdges; % these edges are not optimized

removingActive  = optData.removingActive;
addingActive    = optData.addingActive;
% 0 ~ ALL, 1 ~ all REM, BND ADD, 2 ~ BND
edgesToCheck    = optData.edgesToCheck;

% Generate random seed (randomly generate edges removed from beginning):
nPixels = randi(OP.BF.nUnknowns, 1);
globIndsRem = setdiff(sort(randperm(OP.BF.nUnknowns, nPixels)), protectedEdges);
globIndsBF  = setdiff(1:OP.BF.nUnknowns, globIndsRem);

% It can be done differently, e.g.:
% globIndsBF = round(rand(1, OP.BF.nUnknowns));
% globIndsBF(protectedEdges) = true;
% globIndsBF = find(globIndsBF);

%% Run greedy step from given point for topology significant edges only
t0 = tic;
% Initial structure (admittance matrix)
Y = inv(OP.Zsystem(globIndsBF,globIndsBF));

GenusInit = false(OP.BF.nUnknowns, 1);
GenusInit(globIndsBF) = true;

% Initial step to begin with:
I  = Y * OP.V(globIndsBF);
FF = feval(fitnessFunction, OP, I, globIndsBF, [], 0);

%% ========================================================================
iter     = 1;
tt(iter) = toc(t0);
minDffTopomin = -inf;
candidates = [0 0];
minDffTopoRem = inf; ffTopoRem = [];
minDffTopoAdd = inf; ffTopoAdd = [];
GenusHistory = zeros(1, OP.BF.nUnknowns);
while minDffTopomin < 0
    if edgesToCheck == 0 % Use all edges
        globIndsBFRemTest = setdiff(globIndsBF, protectedEdges);
        globIndsBFAddTest = setdiff(1:OP.BF.nUnknowns, globIndsBF);
    elseif edgesToCheck == 1 % Use only boundary for just ADDING
        [~, globIndsBFAddTest] = objectBoundary(OP.BF, globIndsBF);
        globIndsBFRemTest = setdiff(globIndsBF, protectedEdges);
        globIndsBFAddTest = setdiff(globIndsBFAddTest, protectedEdges);
    elseif edgesToCheck == 2 % Use boundary both for ADDING and REMOVING
        [globIndsBFRemTest, globIndsBFAddTest] = objectBoundary(OP.BF, globIndsBF);
        globIndsBFRemTest = setdiff(globIndsBFRemTest, protectedEdges);
        globIndsBFAddTest = setdiff(globIndsBFAddTest, protectedEdges);
    end
    
    if removingActive
        % Evaluate all removals
        IB_rem = topoRemove(Y, I, globIndsBF, globIndsBFRemTest);
        % Evaluate topology sensitivity (for removal)
        ffTopoRem = feval(fitnessFunction, ...
            OP, IB_rem, globIndsBF, globIndsBFRemTest, -1);
        % Find the most valuable edge to be removed
        [minDffTopoRem, locIndBFrem] = min(ffTopoRem - FF(iter));
    end
    
    if addingActive
        % Evaluate all additionals (USE globIndsBF_YB as output if needee)!
        IB_add = topoAdd(OP.Zsystem, OP.V, Y, I, globIndsBF, globIndsBFAddTest);
        % Evaluate topology sensitivity (for removal)
        ffTopoAdd = feval(fitnessFunction, ...
            OP, IB_add, globIndsBF, globIndsBFAddTest, +1);
        % Find the most valuable edge to be added
        [minDffTopoAdd, locIndBFadd] = min([ffTopoAdd - FF(iter); inf]);
    end
    candidates = candidates + [length(ffTopoRem) length(ffTopoAdd)];
    
    % Evaluate gradients of topology senstivities...
    minDffTopomin = min(minDffTopoRem, minDffTopoAdd);
    if minDffTopomin < 0
        iter = iter + 1;
        FF(iter) = FF(iter-1) + minDffTopomin;
        if minDffTopoRem < minDffTopoAdd % REMOVE EDGE
            % Update current and remove zero item
            locIndIrem = glob2loc(globIndsBF, globIndsBFRemTest(locIndBFrem));
            I = IB_rem([1:(locIndIrem-1), (locIndIrem+1):end], locIndBFrem);
            % Update indices and admittance matrix            
            globIndBFrem = globIndsBFRemTest(locIndBFrem);
            [Y, globIndsBF] = removeEdges(Y, globIndsBF, globIndBFrem);
            % Record to history
            GenusHistory(iter-1) = -globIndBFrem;
        else % ADD EDGE
            % Update current and sort it properly
            [~, inds] =  sort([globIndsBF, globIndsBFAddTest(locIndBFadd)]);
            I = IB_add(inds, locIndBFadd);
            % Update indices and admittance matrix
            globIndBFadd = globIndsBFAddTest(locIndBFadd);
            [Y, globIndsBF] = addEdges(Y, OP.Zsystem, globIndsBF, globIndBFadd);
            % Record to history            
            GenusHistory(iter-1) = +globIndBFadd;
        end
    end
    fprintf(1, 'Iter: %04d, abs improv = %1.5e, relativ. improv = %1.5e, type: %d\n', ...
        iter-1, minDffTopomin, minDffTopomin/FF(iter-1), GenusHistory(iter-1));
    tt(iter) = toc(t0);
end

GenusFinal = false(OP.BF.nUnknowns, 1);
GenusFinal(globIndsBF) = true;