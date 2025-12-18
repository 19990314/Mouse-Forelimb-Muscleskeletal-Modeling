% mouse_reach_Scurve_MC_oneTarget_fromWorkspace.m
%
% Purpose:
%   - Load from mouse_workspace_samples.mat:
%       modelPath, coordNames, markerName, q_start, q_end_tasks, q_min, q_max
%   - Pick one target k_demo:
%       Baseline trajectory = standard S-curve (q_start -> q_end_tasks(k_demo,:))
%   - Monte Carlo:
%       Use makeRandomFamilyTrajectory to generate a set of "family" trajectories
%       For each trajectory -> OpenSim -> compute all muscles' V_norm(t)
%       Compute max |V_norm|, plot histograms, and the proportion of "bad" trajectories
%       ★ New: Plot 3D end-effector trajectories (baseline + MC family),
%              all translated so the start point is the origin

clear; clc; close all;
import org.opensim.modeling.*;

%% 1. Load workspace results
dataFile = 'mouse_workspace_samples.mat';
if ~isfile(dataFile)
    error('File "%s" not found. Please run mouse_workspace_sampling.m first.', dataFile);
end
S = load(dataFile);

modelPath   = S.modelPath;
coordNames  = S.coordNames;
markerName  = S.markerName;   % usually 'handm'
q_start     = S.q_start;
q_end_tasks = S.q_end_tasks;
q_min       = S.q_min;
q_max       = S.q_max;

K  = size(q_end_tasks,1);
nq = numel(coordNames);

fprintf('Loaded workspace from %s, K = %d targets.\n', dataFile, K);

% target
k_demo = 100;
q_end  = q_end_tasks(k_demo,:)';
fprintf('Using target %d as demo (q_end).\n', k_demo);

fprintf('Loading model: %s\n', modelPath);
model  = Model(modelPath);
state  = model.initSystem();

coordSet = model.getCoordinateSet();
muscles  = model.getMuscles();
nm       = muscles.getSize();

lopt        = zeros(1,nm);
muscleNames = cell(nm,1);
for m = 0:nm-1
    mus = muscles.get(m);
    lopt(m+1)        = mus.getOptimalFiberLength();
    muscleNames{m+1} = char(mus.getName());
end

fprintf('Model has %d muscles.\n', nm);

T_total = 0.3;
N_time  = 201;                        % odd number is convenient for finite differences
time    = linspace(0, T_total, N_time)';

% Baseline joint trajectory
[Q_s, Qdot_s] = makeSCurveTrajectory(time, q_start, q_end);
[~, Vm_s]     = computeMuscleKinematicsForTrajectory( ...
                    model, state, coordNames, time, Q_s, Qdot_s);
Vm_s_norm     = Vm_s ./ lopt;

% ★ Baseline end-effector trajectory (absolute coordinates)
xyz_s_abs = computeEndEffTrajectoryFromMarker( ...
    model, state, coordNames, time, Q_s, markerName);

% Translate so baseline start point is the origin
origin    = xyz_s_abs(1,:);                % [1x3]
xyz_s_rel = xyz_s_abs - origin;            % baseline in relative frame

%% 4. Monte Carlo
nMC      = 50;
v_thresh = 5.0;

fprintf('\n==================== MONTE CARLO AROUND S-CURVE ====================\n');
fprintf('Target index = %d\n', k_demo);
fprintf('nMC = %d, v_thresh = %.2f Lopt/s\n\n', nMC, v_thresh);

% Random seed
tstart = clock;
seed   = round(sum(1000*tstart));
rand('state', seed);

badCount      = 0;
runningBad    = [];
goodIndex     = [];
coefvariation = [];

peakSpeed_all = zeros(nMC, nm);
meanSpeed_all = zeros(nMC, nm);

% ★ Store end-effector relative trajectories for each MC trial, for 3D plotting
xyz_MC_rel = cell(nMC,1);

for i = 1:nMC
    fprintf('MC trajectory %d / %d ...\n', i, nMC);

    % 1) Random family trajectory (start/end still q_start, q_end)
    [Q_i, Qdot_i] = makeRandomFamilyTrajectory(time, q_start, q_end);

    % 2) Clamp to ROM
    for j = 1:nq
        Q_i(:,j) = min(Q_i(:,j), q_max(j));
        Q_i(:,j) = max(Q_i(:,j), q_min(j));
    end

    % 3) Muscle kinematics
    [~, Vm_i]  = computeMuscleKinematicsForTrajectory( ...
                     model, state, coordNames, time, Q_i, Qdot_i);
    Vm_i_norm  = Vm_i ./ lopt;

    peakSpeed_all(i,:) = max(abs(Vm_i_norm), [], 1);
    meanSpeed_all(i,:) = mean(abs(Vm_i_norm), 1);

    maxV_this = max(peakSpeed_all(i,:));

    if maxV_this > v_thresh
        badCount   = badCount + 1;
        runningBad = [runningBad; badCount / i];
        goodIndex  = [goodIndex; i];

        coefvariation = [coefvariation; ...
                         100*std(runningBad)/mean(runningBad)];
    end

    % 4) ★ End-effector 3D trajectory for this MC trial (absolute coordinates)
    xyz_i_abs = computeEndEffTrajectoryFromMarker( ...
        model, state, coordNames, time, Q_i, markerName);

    % Translate to the coordinate frame whose origin is the baseline start point
    xyz_MC_rel{i} = xyz_i_abs - origin;
end

if isempty(runningBad)
    fprintf('\nAmong nMC = %d family trajectories, max |V_norm| never exceeded %.2f Lopt/s.\n', ...
        nMC, v_thresh);
else
    fprintf('\nMonte Carlo results:\n');
    fprintf('  Number of bad trajectories = %d / %d (%.2f %%)\n', ...
        badCount, nMC, 100*badCount/nMC);
end

%% 5. Plot runningBad & CV (same as original script)

if ~isempty(runningBad)
    figure('Name','MC proportion of bad trajectories (around S-curve)');
    subplot(2,1,1);
    plot(goodIndex, runningBad, '.','MarkerSize',10);
    xlabel('Trajectory index i');
    ylabel('Proportion of bad trajectories');
    title(sprintf('Progression (v_{th} = %.2f Lopt/s)', v_thresh));
    grid on;

    subplot(2,1,2);
    plot(goodIndex, coefvariation, '.','MarkerSize',10);
    xlabel('Trajectory index i');
    ylabel('CV of proportion (%%)');
    grid on;
end

%% 6. Example: histogram of max/mean V_norm for one muscle

mPlot = 20;  % e.g., set to the index for Biceps or Pec as you like
muscleName = muscleNames{mPlot};

figure('Name','Histogram of Vnorm for one muscle');
subplot(2,1,1);
histogram(peakSpeed_all(:,mPlot));
xlabel('max |V_{norm}| (L_{opt}/s)');
ylabel('Count');
title(sprintf('Muscle %s: max V_{norm}', muscleName), 'Interpreter','none');
grid on;

subplot(2,1,2);
histogram(meanSpeed_all(:,mPlot));
xlabel('mean |V_{norm}| (L_{opt}/s)');
ylabel('Count');
title(sprintf('Muscle %s: mean |V_{norm}|', muscleName), 'Interpreter','none');
grid on;

%% 7. Ranking: which muscles have the highest average peak speed?

meanPeak = mean(peakSpeed_all, 1);
[sortedVal, idxSorted] = sort(meanPeak, 'descend');

fprintf('\nMuscles ranked by mean max |V_norm| over MC family trajectories:\n');
for k = 1:nm
    fprintf('%2d: %-30s  %.3f Lopt/s\n', ...
        k, muscleNames{idxSorted(k)}, sortedVal(k));
end

%% 8. ★ 3D plot: baseline + Monte Carlo family end-effector trajectories (in mm)

nPlot = min(50, nMC);   % plot at most 50 family trajectories to avoid clutter
mm = 1000;              % m -> mm

figure('Name','End-effector 3D trajectories - baseline + MC family (mm)');
hold on; grid on; axis equal;

% baseline S-curve (thick black line)  [mm]
plot3(mm*xyz_s_rel(:,1), mm*xyz_s_rel(:,2), mm*xyz_s_rel(:,3), ...
      'k-', 'LineWidth', 2);

% Monte Carlo family trajectories (thin lines)  [mm]
for i = 1:nPlot
    xyz_i_rel = xyz_MC_rel{i};  % still meters in storage
    plot3(mm*xyz_i_rel(:,1), mm*xyz_i_rel(:,2), mm*xyz_i_rel(:,3), ...
          'LineWidth', 0.8);
end

xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title(sprintf('Baseline S-curve + %d MC family trajectories (start at [0 0 0] mm)', nPlot));
legend({'Baseline S-curve','MC family'}, 'Location','best');
view(45, 20);
hold off;


%% ===== helper functions =====

function [Q, Qdot] = makeSCurveTrajectory(time, q_start, q_end)
N  = numel(time);
nq = numel(q_start);
Q  = zeros(N, nq);

tau = (time - time(1)) / (time(end) - time(1));   % 0 -> 1
s   = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;

for k = 1:N
    Q(k,:) = (q_start + (q_end - q_start) * s(k)).';
end

Qdot = finiteDiff(time, Q);
end

function [Q, Qdot] = makeRandomFamilyTrajectory(time, q_start, q_end)
N  = numel(time);
nq = numel(q_start);
Q  = zeros(N, nq);

T = time(end) - time(1);

for i = 1:nq
    Ti_factor = 0.5 + 0.5*rand();    % [0.5, 1.0]
    Ti        = T * Ti_factor;

    tau = (time - time(1)) / Ti;
    tau(tau > 1) = 1;

    s_base = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;

    alpha = 0.2 + 0.6*rand();
    sigma = 0.1 + 0.1*rand();
    bump  = exp(-((tau - alpha)/sigma).^2);

    eps_amp = 0.1 * (2*rand - 1);
    s_raw   = s_base + eps_amp * bump;

    s_min = s_raw(1);
    s_max = s_raw(end);
    if abs(s_max - s_min) < 1e-6
        s = s_base;
    else
        s = (s_raw - s_min) / (s_max - s_min);
    end

    Q(:,i) = q_start(i) + (q_end(i) - q_start(i)) * s;
end

Qdot = finiteDiff(time, Q);
end

function Qdot = finiteDiff(time, Q)
[N, nq] = size(Q);
Qdot    = zeros(N, nq);

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

function [Lm, Vm] = computeMuscleKinematicsForTrajectory( ...
    model, state, coordNames, time, Q, Qdot)

import org.opensim.modeling.*;

coordSet = model.getCoordinateSet();
muscles  = model.getMuscles();

N  = size(Q,1);
nq = numel(coordNames);
nm = muscles.getSize();

Lm = zeros(N, nm);
Vm = zeros(N, nm);

for k = 1:N
    for j = 1:nq
        coord = coordSet.get(coordNames{j});
        coord.setValue(state,      Q(k,j));
        coord.setSpeedValue(state, Qdot(k,j));
    end

    model.realizeVelocity(state);

    for m = 0:nm-1
        musc = muscles.get(m);
        Lm(k, m+1) = musc.getLength(state);
        Vm(k, m+1) = musc.getFiberVelocity(state);
    end
end
end

function xyz = computeEndEffTrajectoryFromMarker( ...
    model, state, coordNames, time, Q, markerName)
% Compute the marker's 3D trajectory in ground coordinates

import org.opensim.modeling.*;

coordSet  = model.getCoordinateSet();
markerSet = model.getMarkerSet();
engine    = model.getSimbodyEngine();

marker      = markerSet.get(markerName);
parentFrame = marker.getParentFrame();
locInParent = marker.get_location();   % Vec3

N  = size(Q,1);
nq = numel(coordNames);
xyz = zeros(N,3);

for k = 1:N
    for j = 1:nq
        coord = coordSet.get(coordNames{j});
        coord.setValue(state, Q(k,j));
    end

    model.realizePosition(state);

    pGround = Vec3();
    engine.getPosition(state, parentFrame, locInParent, pGround);

    xyz(k,1) = pGround.get(0);
    xyz(k,2) = pGround.get(1);
    xyz(k,3) = pGround.get(2);
end
end
