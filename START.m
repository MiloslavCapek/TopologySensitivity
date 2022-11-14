%% Script performing local step of topology sensitivity based on exact 
%  reanalysis
% 
% The script starts from a random seed (Monte Carlo) and repeats local updates
% as long as the local minimum is reached (i.e., no addition or removal can
% minimize the fitness function further).
% 
% The additions and removals are using Sherman-Morrison-Woodbury identity 
% and are fully vectorized.
% 
% This procedure is general and can be used to optimize any antenna metric 
% (or combination of more metrics). In its full version, it is combined 
% with heuristic (global) step maintaining diversity and robustness, 
% see [1]. Some applications are in [2].
%
% The method uses matrices stemming from method-of-moments solution to field 
% integral equation. In this script, the precalculated matrices results from 
% electric field integral equation. Rao-Wilton-Glisson basis functions
% are used as optimization degrees of freedom.
% 
% [1] Capek, M., Jelinek, L., Kadlec, P., Gustafsson, M.: Memetic Scheme 
%     for Inverse Design Using an Exact Reanalysis of Method-of-Moments
%     Models - Part 1: Theory and Implementation, submitted, 
%     pp. 1-12, 2022.
% [2] Capek, M., Jelinek, L., Kadlec, P., Gustafsson, M.: Memetic Scheme 
%     for Inverse Design Using an Exact Reanalysis of Method-of-Moments
%     Models - Part 2: Examples and Properties, pp. 1-13, 2022.
% 
% 2022, Miloslav Capek, CTU in Prague, miloslav.capek@fel.cvut.cz

clear;
clc;

% Precalculated matrices are loaded (they can be changed in AToM [1]):
load('plate_OP_structure.mat');

% Auxiliarly matrices
OP.omXm    = omXm; % stored magnetic energy matrix
OP.omXe    = omXe; % stored electric energy matrix
OP.Zsystem = Z;    % matrix used for rank-1 additions/removals and updates
OP.Z       = Z;    % impedance matrix used to evaluated radiated power
OP.S       = S1;   % real(Z) = S1'*S1 decomposition of radiation matrix
OP.V       = V;    % excitation vector
OP.Mesh    = Mesh; % MATLAB structure of geometry (from AToM)
OP.BF      = BF;   % MATLAB structure of basis functions (from AToM)

%% Fitness function settings (excitation, type of problem)
OP.ports        = port;

% Q-factor minimization without penalty imposing self-resonance
optData.fitness = @ff_minQ;

% Q-factor minimization with penalty imposing self-resonance
% optData.fitness = @ff_minQ_selfres;

% Optimization settings:
optData.protectedEdges = OP.ports; % these edges are not optimized
optData.removingActive = true;     % can be edges removed (?)
optData.addingActive   = true;     % can be edges adeed (?)
optData.edgesToCheck   = 0;        % 0 ~ ALL, 1 ~ all REM, BND ADD, 2 ~ BND

%% Performs series of local steps to reach the nearest local minimum (from 
% a random seed)
t0 = tic;
% Here, the optimization starts from a random seed:
[FF, GenusInit, GenusFinal, GenusHistory, tt, candidates] = ...
    localSteps(OP, optData);
toc(t0)

% Results:
Qmin     = FF(end)
QminNorm = Qmin * ka^3 % here, Q-factor is normalized to electrical size
% To evaluate lower bound to Q-factor, use package FunBo: 
% http://antennatoolbox.com/fundamentalBounds

%% Time course of the cost function starting from a random seed and finishing
% at local minimum
figure;
% loglog(FF, '.-'); xlabel('local step iteration (-)');
loglog(tt, FF, 'x-'); xlabel('computational time (t)');
grid on;
ylabel('Q (-)');

%% Show connectivity of the locally optimal design
figure('color', 'w', 'pos', [250 250 800 600]);
trisurf(OP.Mesh.connectivityList, ...
    OP.Mesh.nodes(:,1), OP.Mesh.nodes(:,2), OP.Mesh.nodes(:,3), ...
    'FaceColor', 3/4*[1 1 1]);
view(0, 90);
axis equal;
for n = 1:OP.BF.nUnknowns
    Ppt = OP.Mesh.triangleCentroids(OP.BF.data(n, 2), :);
    Mpt = OP.Mesh.triangleCentroids(OP.BF.data(n, 4), :);
    if GenusFinal(n)
        if n == OP.ports
            opt = {'LineWidth', 10, 'Color', 'b'};
        else
            opt = {'LineWidth', 3, 'Color', 'r'};
        end
        line('XData', [Ppt(1) Mpt(1)], ...
             'YData', [Ppt(2) Mpt(2)], ...
             'ZData', [Ppt(3) Mpt(3)], opt{:});
    end
end
xlim([min(OP.Mesh.nodes(:,1)) max(OP.Mesh.nodes(:,1))]);
ylim([min(OP.Mesh.nodes(:,2)) max(OP.Mesh.nodes(:,2))]);

%% Show locally optimal current density (requires AToM)
I = zeros(OP.BF.nUnknowns, 1);
I(GenusFinal) = OP.Zsystem(GenusFinal, GenusFinal) \ OP.V(GenusFinal);
try
    hndl = results.plotCurrent(OP.Mesh, 'basisFcn', OP.BF, ...
        'Ivec', I, 'part', 'abs', ...
        'arrowScale', 'proportional');
catch
    fprintf(1, ['You have to install AToM first to see the current.\n', ...
            'Visit: <a href = "http://antennatoolbox.com">', ...
            'antennatoolbox.com</a> and download FREE version.\n']);
end
