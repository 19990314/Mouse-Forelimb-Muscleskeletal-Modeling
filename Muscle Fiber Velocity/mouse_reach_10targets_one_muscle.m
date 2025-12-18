clear; clc; close all;
import org.opensim.modeling.*;

%% ------------------------------------------------------------------------
% 1. Load workspace
% -------------------------------------------------------------------------
dataFile = 'mouse_workspace_samples.mat';
if ~isfile(dataFile)
    error('File "%s" not found. RUN mouse_workspace_sampling.m', dataFile);
end

S = load(dataFile);

modelPath     = S.modelPath;
coordNames    = S.coordNames;
markerName    = S.markerName;
q_start       = S.q_start;          % [6 x 1] rad
q_end_tasks   = S.q_end_tasks;      % [K x 6] rad
xyz_end_tasks = S.xyz_end_tasks;    % [K x 3]
startXYZ      = S.startXYZ;         % [1 x 3]
q_min         = S.q_min;            
q_max         = S.q_max;            

K  = size(q_end_tasks,1);
nq = numel(coordNames);

fprintf('Loaded workspace from %s\n', dataFile);
fprintf('  Number of targets K = %d\n', K);
fprintf('  Using marker: %s\n\n', markerName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Target Switch %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
onlyInspect = false;     % <<< true = only selected target；false = all targets
kInspect    = 100;      % <<< tar target（1~K）

if onlyInspect
    if kInspect < 1 || kInspect > K
        error('kInspect=%d out of range. K=%d', kInspect, K);
    end
    kList = kInspect;          % only selected one Ex. #100
else
    kList = 1:K;               % All
end
K_run = numel(kList);          

%% ------------------------------------------------------------------------
% 2. Time line
% -------------------------------------------------------------------------
T_total = 0.3;                  % s
N_time  = 201;
time    = linspace(0, T_total, N_time)';   % [N_time x 1]

%% ------------------------------------------------------------------------
% 3. Load model + muscles
% -------------------------------------------------------------------------
fprintf('Loading model: %s\n', modelPath);
model = Model(modelPath);
state = model.initSystem();

coordSet = model.getCoordinateSet();
muscles  = model.getMuscles();
nm       = muscles.getSize();
fprintf('Model has %d muscles.\n', nm);

lopt        = zeros(1, nm);
muscleNames = cell(nm,1);
for m = 0:nm-1
    musc = muscles.get(m);
    lopt(m+1)        = musc.getOptimalFiberLength();
    muscleNames{m+1} = char(musc.getName());
end

% Color the muscle
cmapBase = [
    0.894 0.102 0.110
    0.215 0.494 0.721
    0.302 0.686 0.290
    0.596 0.306 0.639
    1.000 0.498 0.000
    0.651 0.337 0.157
    0.969 0.506 0.749
    0.121 0.470 0.705
    0.200 0.627 0.173
    0.890 0.466 0.760
    0.737 0.741 0.133
    0.090 0.745 0.811
    0.556 0.184 0.419
    0.984 0.603 0.600
    0.682 0.780 0.909
    0.992 0.749 0.435
    0.600 0.600 0.600
    0.145 0.203 0.580
    0.305 0.654 0.698
    0.800 0.470 0.000
    0.000 0.000 0.000
];
if nm <= size(cmapBase,1)
    cmapMuscles = cmapBase(1:nm,:);
else
    rep = ceil(nm / size(cmapBase,1));
    cmapMuscles = repmat(cmapBase, rep, 1);
    cmapMuscles = cmapMuscles(1:nm,:);
end

%% ------------------------------------------------------------------------
% 4. Output folder
% -------------------------------------------------------------------------
if onlyInspect
    outDir = sprintf('Scurve_results_Target%03d_only', kInspect); % ★ Only selected
else
    outDir = sprintf('Scurve_results_%d_targets', K);             % ★ All
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
fprintf('All figures and data will be saved in folder: %s\n\n', outDir);

%% ------------------------------------------------------------------------
% 5. Compute S-curve + muscle kinematics (only kList)
% -------------------------------------------------------------------------
N = N_time;

xyz_all_abs       = cell(K_run,1);
maxV_by_target    = zeros(K_run, nm);          % [K_run x nm]
Lm_all_targets    = zeros(N_time, nm, K_run); % [N_time x nm x K_run]
dL_net_by_target  = zeros(K_run, nm);         % [K_run x nm]

fprintf('Computing S-curve trajectories and muscle kinematics for %d target(s)...\n', K_run);

for ii = 1:K_run
    k = kList(ii);                  %  target selected（example: 100）
    q_end = q_end_tasks(k,:)';       % [6 x 1] rad

    % ---- 5.1 S-curve trajectory
    [Q_k, Qdot_k] = makeSCurveTrajectory(time, q_start, q_end);

    % ---- 5.2 muscle kinematics
    [Lm_k, Vm_k] = computeMuscleKinematicsForTrajectory( ...
        model, state, coordNames, time, Q_k, Qdot_k);

    Vm_norm_k = Vm_k ./ lopt;        % [N_time x nm] (Lopt/s)

    %
    Lm_all_targets(:,:,ii) = Lm_k;

    dL_net = Lm_k(end,:) - Lm_k(1,:);            % [1 x nm]
    dL_net_by_target(ii,:) = dL_net;

    maxV_by_target(ii,:) = max(abs(Vm_norm_k), [], 1);

    % ---- 5.3 end-effector (marker) trajectory in ground
    xyz_abs_k = computeEndEffTrajectoryFromMarker( ...
        model, state, coordNames, Q_k, markerName);
    xyz_all_abs{ii} = xyz_abs_k;

    % ---- 5.4 Vnorm(t): all muscles
    figTarget = figure('Name', sprintf('All muscles V_{norm} - Target %d', k), ...
                       'Visible','off');
    hold on;
    for m = 1:nm
        plot(time, Vm_norm_k(:,m), ...
             'Color', cmapMuscles(m,:), ...
             'LineWidth', 1.3);
    end
    yline(0,'k--');
    xlabel('Time (s)');
    ylabel('V_{norm} (L_{opt}/s)');
    title(sprintf('All muscles normalized velocity - Target %d', k));
    grid on;
    legend(muscleNames, 'Interpreter','none', 'Location','eastoutside');
    hold off;

    fnamePng = fullfile(outDir, sprintf('Vnorm_allMuscles_Target%03d.png', k));
    saveas(figTarget, fnamePng);
    close(figTarget);

    % ---- 5.5 Top5 |ΔL_net| 的 ΔL(t) + bar
    [~, idxSort] = sort(abs(dL_net), 'descend');
    nTop = min(5, nm);
    idxTop = idxSort(1:nTop);

   % 5.5a Top5 ΔL(t)  (mm)
figTopDL = figure('Name', sprintf('Top%d muscles \\DeltaL(t) - Target %d', nTop, k), ...
                  'Position',[200 200 800 600], ...
                  'Visible','off');
hold on;
for jj = 1:nTop
    mIdx = idxTop(jj);
    dL_t_mm = 1000 * (Lm_k(:,mIdx) - Lm_k(1,mIdx));   % <<< m -> mm
    plot(time, dL_t_mm, ...
         'Color', cmapMuscles(mIdx,:), ...
         'LineWidth', 1.8, ...
         'DisplayName', muscleNames{mIdx});
end
yline(0,'k--');
xlabel('Time (s)');
ylabel('\Delta L (mm)');
title(sprintf('Top %d muscles length change \\DeltaL(t) - Target %d', nTop, k), ...
      'Interpreter','tex');
grid on;
legend('Location','eastoutside','Interpreter','none');

ax = gca;
ax.YAxis.Exponent = 0;   % <<< 不要 ×10^n
hold off;

fnameTopDL = fullfile(outDir, sprintf('Top%d_dL_time_Target%03d_mm.png', nTop, k));
saveas(figTopDL, fnameTopDL);
close(figTopDL);



% 5.5b Net ΔL bar (Top5)
figBar = figure('Name', sprintf('Net \\DeltaL for all muscles - Target %d', k), ...
                'Position',[200 200 900 500], ...
                'Visible','off');
hold on;

dL_net_mm = 1000 * dL_net;   % m -> mm

bar(1:nm, dL_net_mm, 'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none');
bar(idxTop, dL_net_mm(idxTop), 'FaceColor',[0.85 0.1 0.1], 'EdgeColor','none');

xlabel('Muscle');
ylabel('\Delta L_{net} (mm)');
title(sprintf('Net length change \\DeltaL_{net} = L_{end}-L_{start} - Target %d', k), ...
      'Interpreter','tex');
grid on;

set(gca,'XTick',1:nm, ...
        'XTickLabel',muscleNames, ...
        'XTickLabelRotation',45, ...
        'TickLabelInterpreter','none');   % 

ax = gca;
ax.YAxis.Exponent = 0;  % 
hold off;

fnameBar = fullfile(outDir, sprintf('Net_dL_bar_allMuscles_Target%03d_mm.png', k));
saveas(figBar, fnameBar);
close(figBar);


    fprintf('  Saved Target %d figures.\n', k);
end

fprintf('Done per-target figures.\n\n');



%% ------------------------------------------------------------------------
% 6. 3D trajectories (relative to common startXYZ = 0, in mm, NO cloud)
% -------------------------------------------------------------------------
scale = 1000;

fig3D = figure('Name','handm S-curve trajectories (relative)', ...
               'Visible','off', 'Color','w');
ax = axes('Parent',fig3D);
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');

% common origin = startXYZ
plot3(ax, 0,0,0, 'kp','MarkerFaceColor','y','MarkerSize',12);
text(ax, 0,0,0, '  start', 'FontSize',10,'Color','k');

for ii = 1:K_run
    k = kList(ii);

    xyz_abs_k = xyz_all_abs{ii};                 % m
    xyz_plot  = (xyz_abs_k - startXYZ) * scale;  % <<< 统一归0 + mm

    col = [0.75 0.75 0.75]; lw = 1.2;
    if k == 100
        col = [0.85 0.10 0.10]; lw = 3.0;
    end

    plot3(ax, xyz_plot(:,1), xyz_plot(:,2), xyz_plot(:,3), ...
          '-', 'LineWidth', lw, 'Color', col);

    if k == 100
        pEnd = xyz_plot(end,:);
        plot3(ax, pEnd(1), pEnd(2), pEnd(3), 'o', ...
              'MarkerFaceColor','r', 'MarkerSize',8, 'MarkerEdgeColor','k');
        text(ax, pEnd(1), pEnd(2), pEnd(3), sprintf('  Target %d', k), ...
             'FontSize',10, 'Color','k');
    end
end

xlabel(ax,'\Delta X (mm)');
ylabel(ax,'\Delta Y (mm)');
zlabel(ax,'\Delta Z (mm)');
title(ax, sprintf('handm S-curve trajectories (relative to start, %d targets)', K_run), ...
      'Interpreter','none');
view(ax, 45, 25);
hold(ax,'off');

saveas(fig3D, fullfile(outDir, sprintf('handm_Scurve_3D_%dTargets_relativeToStart_mm.png', K_run)));
close(fig3D);



%% ------------------------------------------------------------------------
% 7. Save MAT results for polar / later analysis
% -------------------------------------------------------------------------
resultsMat = fullfile(outDir, 'Scurve_results_for_polar.mat');
save(resultsMat, ...
     'kList', 'onlyInspect', ...
     'time', 'modelPath', 'coordNames', 'markerName', ...
     'q_start', 'q_end_tasks', 'xyz_end_tasks', 'startXYZ', ...
     'muscleNames', 'lopt', ...
     'maxV_by_target', ...
     'Lm_all_targets', 'dL_net_by_target');

fprintf('Numerical results saved to %s\n', resultsMat);
fprintf('All done.\n');

%% ========================================================================
% Local functions
% ========================================================================
function [Q, Qdot] = makeSCurveTrajectory(time, q_start, q_end)
N  = numel(time);
nq = numel(q_start);
Q  = zeros(N, nq);

tau = (time - time(1)) / (time(end) - time(1));   % 0->1
s   = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;

for k = 1:N
    Q(k,:) = (q_start + (q_end - q_start) * s(k)).';
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

function [Lm, Vm] = computeMuscleKinematicsForTrajectory(model, state, coordNames, time, Q, Qdot)
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

function xyz = computeEndEffTrajectoryFromMarker(model, state, coordNames, Q, markerName)
import org.opensim.modeling.*;

coordSet  = model.getCoordinateSet();
markerSet = model.getMarkerSet();
engine    = model.getSimbodyEngine();

marker      = markerSet.get(markerName);
parentFrame = marker.getParentFrame();
locInParent = marker.get_location();

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
