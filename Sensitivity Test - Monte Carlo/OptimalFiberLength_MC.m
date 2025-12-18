% ========================================================================
%  mouse_LoptMonteCarlo_Vnorm.m
% ========================================================================

clear; clc; close all;
import org.opensim.modeling.*;

%% ------------------------------------------------------------------------
% 1. Load workspace info
% -------------------------------------------------------------------------
dataFile = 'mouse_workspace_samples.mat';
if ~isfile(dataFile)
    error('File "%s" not found. Please run mouse_workspace_sampling.m first.', dataFile);
end

S = load(dataFile);
modelPath     = S.modelPath;
coordNames    = S.coordNames;      % names of the 6 DOFs
q_start       = S.q_start;         % [6 x 1] rad
q_end_tasks   = S.q_end_tasks;     % [K x 6] rad
markerName    = S.markerName;      % not used here, but printed for reference

K  = size(q_end_tasks,1);
nq = numel(coordNames);

fprintf('Loaded workspace info from %s\n', dataFile);
fprintf('  Model:   %s\n', modelPath);
fprintf('  Marker:  %s\n', markerName);
fprintf('  DOFs:    %d\n', nq);
fprintf('  Targets: %d\n\n', K);

% Choose one target posture as the S-curve end posture (you can change k_demo)
k_demo = 1;
q_end  = q_end_tasks(k_demo,:)';
fprintf('Using target %d as S-curve end posture.\n\n', k_demo);

%% ------------------------------------------------------------------------
% 2. Time vector + S-curve trajectory (all Lopt variants reuse the same q(t))
% -------------------------------------------------------------------------
T_total = 0.3;          % total duration (s); try 0.3 / 0.8 to see effects on speed
N_time  = 201;          % odd number of points for central differencing
time    = linspace(0, T_total, N_time)';   % [N x 1]

[Q_baseline, Qdot_baseline] = makeSCurveTrajectory(time, q_start, q_end);

%% ------------------------------------------------------------------------
% 3. Baseline model: Lopt_base + baseline max |V_norm|
% -------------------------------------------------------------------------
fprintf('Loading baseline model and computing baseline V_{norm}...\n');

baseModel = Model(modelPath);
baseMuscles = baseModel.getMuscles();
nm = baseMuscles.getSize();

muscleNames = cell(nm,1);
lopt_base   = zeros(1,nm);
for m = 0:nm-1
    mus = baseMuscles.get(m);
    muscleNames{m+1} = char(mus.getName());
    lopt_base(m+1)   = mus.getOptimalFiberLength();
end
fprintf('Model has %d muscles.\n\n', nm);

[Vm_norm_baseline, Vpeak_baseline] = ...
    compute_Vpeak_for_model(baseModel, coordNames, time, ...
                            Q_baseline, Qdot_baseline, lopt_base);

%% ------------------------------------------------------------------------
% 4. Lopt Monte Carlo parameters (BME404 style)
% -------------------------------------------------------------------------
nLopt   = 100;    % number of Lopt samples; try 50 / 200
rangeL  = 0.10;   % Lopt scaling range: ±10%, i.e., [0.9, 1.1]
v_thresh = 5.0;   % "bad parameter set" threshold for V_norm (Lopt/s)

fprintf('=========== Lopt Monte Carlo (no geometry changes) ===========\n');
fprintf('Number of Lopt variants = %d\n', nLopt);
fprintf('Lopt scale factor ~ U[%.2f, %.2f]\n', 1-rangeL, 1+rangeL);
fprintf('Bad-set threshold v_{th} = %.2f Lopt/s\n\n', v_thresh);

Vpeak_all = zeros(nLopt, nm);   % [variant x muscle]

% BME404-style tracking variables
badCount      = 0;      % number of "bad Lopt sets"
runningBad    = [];     % running estimate of P(bad)
goodIndex     = [];     % sample indices where a bad set occurred
coefvariation = [];     % CV(%) of runningBad

% Random seed (reproducibility)
tstart = clock;
seed   = round(sum(1000*tstart));
rng(seed,'twister');
fprintf('Random seed = %d\n\n', seed);

%% ------------------------------------------------------------------------
% 5. Monte Carlo loop: randomize Lopt, then compute V_norm along the same S-curve
% -------------------------------------------------------------------------
for r = 1:nLopt
    fprintf('Variant %d / %d ...\n', r, nLopt);

    % 5.1 Recreate a model from file each time (avoid accumulating changes)
    model   = Model(modelPath);
    muscles = model.getMuscles();

    % 5.2 Randomly assign an Lopt scale factor to each muscle
    lopt_var = zeros(1,nm);
    for m = 0:nm-1
        mus = muscles.get(m);

        % scaleL ~ U[1-rangeL, 1+rangeL]
        scaleL = (1 - rangeL) + (2*rangeL)*rand();
        L0_new = lopt_base(m+1) * scaleL;

        mus.setOptimalFiberLength(L0_new);
        lopt_var(m+1) = L0_new;
    end

    % 5.3 Under this Lopt set, compute max |V_norm| along the same S-curve
    [~, Vpeak_var] = compute_Vpeak_for_model( ...
        model, coordNames, time, Q_baseline, Qdot_baseline, lopt_var);

    Vpeak_all(r,:) = Vpeak_var;   % [1 x nm]

    maxV_this = max(Vpeak_var);

    if maxV_this > v_thresh
        badCount = badCount + 1;
        runningBad(end+1,1) = badCount / r;   %#ok<SAGROW>
        goodIndex(end+1,1)  = r;              %#ok<SAGROW>

        if numel(runningBad) > 1
            coefvariation(end+1,1) = 100 * std(runningBad) / mean(runningBad); %#ok<SAGROW>
        else
            coefvariation(end+1,1) = NaN;  % cannot compute CV for a single point
        end
    end
end

fprintf('\nAll Lopt variants done.\n\n');

%% ------------------------------------------------------------------------
% 6. BME404-style results: estimated probability of "bad" Lopt sets + convergence
% -------------------------------------------------------------------------
if isempty(runningBad)
    fprintf('Note: among %d Lopt samples, max |V_{norm}| never exceeded %.2f Lopt/s.\n', ...
        nLopt, v_thresh);
else
    finalProb = runningBad(end);
    fprintf('Estimated probability of "bad" Lopt sets ≈ %.2f %% (threshold = %.2f Lopt/s)\n\n', ...
        100*finalProb, v_thresh);

    figure('Name','Lopt Monte Carlo - probability of bad L_{opt} sets');
    subplot(2,1,1);
    plot(goodIndex, runningBad, '.-');
    xlabel('Sample index (r)');
    ylabel('P(bad L_{opt})');
    title('Running estimate of probability of bad L_{opt} sets');
    grid on;

    subplot(2,1,2);
    plot(goodIndex, coefvariation, '.-');
    xlabel('Sample index (r)');
    ylabel('CV of running estimate (%)');
    title('Coefficient of variation of P(bad L_{opt})');
    grid on;
end

%% ------------------------------------------------------------------------
% 7. Sensitivity per muscle: CV of max |V_norm| across Lopt variants
% -------------------------------------------------------------------------
Vmean = mean(Vpeak_all, 1);      % [1 x nm]
Vstd  = std(Vpeak_all, 0, 1);    % [1 x nm]

CV_V = zeros(1,nm);
for m = 1:nm
    if Vmean(m) ~= 0
        CV_V(m) = 100 * Vstd(m) / abs(Vmean(m));
    else
        CV_V(m) = 0;
    end
end

[CV_V_sorted, idxV] = sort(CV_V, 'descend');

fprintf('===== Muscles ranked by CV of max |V_{norm}| across Lopt variants =====\n');
nShow = min(10, nm);
for k = 1:nShow
    m = idxV(k);
    fprintf('%2d: %-30s  CV_V = %6.2f %%   (mean peak = %.3f Lopt/s, baseline = %.3f)\n', ...
        k, muscleNames{m}, CV_V_sorted(k), Vmean(m), Vpeak_baseline(m));
end
fprintf('\n');

%% ------------------------------------------------------------------------
% 8. Bar plot: sensitivity (Top 10 muscles)
% -------------------------------------------------------------------------
figure('Name','Top 10 muscles by V_{norm} sensitivity to L_{opt}');
bar(CV_V(idxV(1:nShow)));
set(gca, 'XTick', 1:nShow, ...
         'XTickLabel', muscleNames(idxV(1:nShow)), ...
         'XTickLabelRotation', 45);
ylabel('CV of max |V_{norm}| (%)');
title(sprintf('Sensitivity of V_{norm} to L_{opt} (\\pm%.0f%% perturbation)', rangeL*100));
grid on;

%% ------------------------------------------------------------------------
% 9. Save results to .mat for later analysis/plotting
% -------------------------------------------------------------------------
LoptMC = struct();
LoptMC.nLopt        = nLopt;
LoptMC.rangeL       = rangeL;
LoptMC.v_thresh     = v_thresh;
LoptMC.Vpeak_all    = Vpeak_all;
LoptMC.Vpeak_baseline = Vpeak_baseline;
LoptMC.Vmean        = Vmean;
LoptMC.CV_V         = CV_V;
LoptMC.badCount     = badCount;
LoptMC.runningBad   = runningBad;
LoptMC.goodIndex    = goodIndex;
LoptMC.coefvariation = coefvariation;
LoptMC.muscleNames  = muscleNames;
LoptMC.coordNames   = coordNames;
LoptMC.time         = time;

save('mouse_LoptMonteCarlo_results.mat', ...
    'modelPath', 'coordNames', 'muscleNames', 'lopt_base', ...
    'Q_baseline', 'Qdot_baseline', 'Vm_norm_baseline', ...
    'LoptMC');

fprintf('All results saved to mouse_LoptMonteCarlo_results.mat\n');
fprintf('Done.\n');

%% ========================================================================
% Helper functions
% ========================================================================

function [Vm_norm, Vpeak] = compute_Vpeak_for_model( ...
    model, coordNames, time, Q_traj, Qdot_traj, lopt)

% Given a model (with Lopt already set) and a trajectory Q(t),
% compute all muscles' normalized fiber velocity:
%   Vm_norm(t) = V_fiber / Lopt
% and
%   Vpeak(m) = max_t |Vm_norm(t,m)|

import org.opensim.modeling.*;

state = model.initSystem();

[~, Vm] = computeMuscleKinematicsForTrajectory( ...
    model, state, coordNames, time, Q_traj, Qdot_traj);

Vm_norm = Vm ./ lopt;                % [N_time x nm]
Vpeak   = max(abs(Vm_norm), [], 1);  % [1 x nm]
end

% -------------------------------------------------------------------------

function [Q, Qdot] = makeSCurveTrajectory(time, q_start, q_end)
% Joint-space quintic S-curve (all DOFs share the same s(t))

N  = numel(time);
nq = numel(q_start);
Q  = zeros(N, nq);

tau = (time - time(1)) / (time(end) - time(1));   % 0 -> 1
s   = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;           % quintic S-curve

for k = 1:N
    Q(k,:) = (q_start + (q_end - q_start) * s(k)).';
end

Qdot = finiteDiff(time, Q);
end

% -------------------------------------------------------------------------

function Qdot = finiteDiff(time, Q)
% Simple central-difference estimate of joint angular velocity

[N, nq] = size(Q);
Qdot = zeros(N, nq);

for j = 1:nq
    qj = Q(:,j);
    for k = 2:N-1
        dt = time(k+1) - time(k-1);
        Qdot(k,j) = (qj(k+1) - qj(k-1)) / dt;
    end
    Qdot(1,j) = (qj(2)   - qj(1))   / (time(2)   - time(1));
    Qdot(N,j) = (qj(N)   - qj(N-1)) / (time(N)   - time(N-1));
end
end

% -------------------------------------------------------------------------

function [Lm, Vm] = computeMuscleKinematicsForTrajectory( ...
    model, state, coordNames, time, Q, Qdot)
% For a trajectory Q(t), compute all muscles' Lm(t) and Vm(t)

import org.opensim.modeling.*;

coordSet = model.getCoordinateSet();
muscles  = model.getMuscles();

N  = size(Q,1);
nq = numel(coordNames);
nm = muscles.getSize();

Lm = zeros(N, nm);
Vm = zeros(N, nm);

for k = 1:N
    % Set joint angles and speeds for this frame
    for j = 1:nq
        coord = coordSet.get(coordNames{j});
        coord.setValue(state,      Q(k,j));
        coord.setSpeedValue(state, Qdot(k,j));
    end

    model.realizeVelocity(state);

    for m = 0:nm-1
        mus = muscles.get(m);
        Lm(k, m+1) = mus.getLength(state);
        Vm(k, m+1) = mus.getFiberVelocity(state);
    end
end
end
