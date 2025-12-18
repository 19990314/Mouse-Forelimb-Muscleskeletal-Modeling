clear; clc; close all;
import org.opensim.modeling.*;

%% 1. Load Workspace
dataFile = 'mouse_workspace_samples.mat';
if ~isfile(dataFile)
    error('File "%s" not found. RUN mouse_workspace_sampling.m', dataFile);
end

S = load(dataFile);

modelPath     = S.modelPath;
coordNames    = S.coordNames;
markerName    = S.markerName;       %'handm'
q_start       = S.q_start;          % [6x1] rad
q_end_tasks   = S.q_end_tasks;      % [K x 6] rad
xyz_end_tasks = S.xyz_end_tasks;    % [K x 3]
startXYZ      = S.startXYZ;         % [1x3]
q_min         = S.q_min;            
q_max         = S.q_max;           

K  = size(q_end_tasks,1);

fprintf('Loaded workspace from %s\n', dataFile);
fprintf('  Number of targets K = %d\n', K);

%% 2. Timeline
T_total = 0.3;          % 总时长 (s)
N_time  = 201;          % 时间采样点数（奇数方便有限差分）
time    = linspace(0, T_total, N_time)';   % [N x 1]

%% 3. Load Muscle parameters
fprintf('Loading model: %s\n', modelPath);
model = Model(modelPath);
state = model.initSystem();

coordSet  = model.getCoordinateSet();
muscles   = model.getMuscles();
nm        = muscles.getSize();

fprintf('Model has %d muscles.\n', nm);

lopt        = zeros(1, nm);
muscleNames = cell(nm,1);
for m = 0:nm-1
    musc = muscles.get(m);
    lopt(m+1)        = musc.getOptimalFiberLength();
    muscleNames{m+1} = char(musc.getName());
end
fprintf('Collected optimal fiber lengths for all muscles.\n\n');

% Color
cmapBase = [
    0.1216 0.4667 0.7059
    0.6824 0.7804 0.9098
    1.0000 0.4980 0.0549
    1.0000 0.7333 0.4706
    0.1725 0.6275 0.1725
    0.5961 0.8745 0.5412
    0.8392 0.1529 0.1569
    1.0000 0.5961 0.5882
    0.5804 0.4039 0.7412
    0.7725 0.6902 0.8353
    0.5490 0.3373 0.2941
    0.7686 0.6118 0.5804
    0.8902 0.4667 0.7608
    0.9686 0.7137 0.8235
    0.4980 0.4980 0.4980
    0.7804 0.7804 0.7804
    0.7373 0.7412 0.1333
    0.8588 0.8588 0.5529
    0.0902 0.7451 0.8118
    0.6196 0.8549 0.8980
    0.8588 0.3490 0.5804
];
if nm <= size(cmapBase,1)
    cmapMuscles = cmapBase(1:nm,:);
else
    rep = ceil(nm / size(cmapBase,1));
    cmapMuscles = repmat(cmapBase, rep, 1);
    cmapMuscles = cmapMuscles(1:nm,:);
end

%% 4.S-curve and calculate normalzied velocity
N = N_time;
Vm_all_norm  = zeros(N, nm, K);  % V_norm(t) [N x nm x K]

fprintf('Computing S-curve trajectories and muscle kinematics for %d targets...\n', K);

for k = 1:K
    q_end_k = q_end_tasks(k,:)';   % [6x1]

    [Q_k, Qdot_k] = makeSCurveTrajectory(time, q_start, q_end_k);

    [~, Vm_k] = computeMuscleKinematicsForTrajectory( ...
        model, state, coordNames, time, Q_k, Qdot_k);

    Vm_k_norm = Vm_k ./ lopt;   % Lopt/s
    Vm_all_norm(:,:,k) = Vm_k_norm;
end
fprintf('Done.\n\n');

%% 6. Arrange the graph
maxVnorm_all = squeeze(max(abs(Vm_all_norm), [], 1));  % [nm x K]
globalMax = max(maxVnorm_all, [], 2);                  % [nm x 1]
[~, idxSort] = sort(globalMax, 'descend');

Ntop = min(8, nm);
idxTop = idxSort(1:Ntop);

% --- rule
isPec = false(size(idxTop));
isTri = false(size(idxTop));
isDel = false(size(idxTop));

for ii = 1:numel(idxTop)
    name = muscleNames{idxTop(ii)};
    isPec(ii) = startsWith(name, 'Pectoralis_');
    isTri(ii) = contains(name, 'Triceps');      % or strcmp(name,'Triceps_Long_Head')
    isDel(ii) = startsWith(name, 'Deltoid_');
end

idxGroup1 = idxTop(isPec | isTri);   % Expectation=4（3 pec + 1 tri）
idxGroup2 = idxTop(~(isPec | isTri));% 

% Check
if numel(idxGroup1) ~= 4 || numel(idxGroup2) ~= 4
    warning('Top8 does not split into 4+4 as expected. Using idxTop(1:4) and idxTop(5:8) fallback.');
    idxGroup1 = idxTop(1:min(4,numel(idxTop)));
    idxGroup2 = idxTop(min(5,numel(idxTop)+1):end);
end

fprintf('Top8 muscles (by global max |Vnorm|):\n');
for ii = 1:Ntop
    m = idxTop(ii);
    fprintf('  #%d: %-30s  globalMax=%.4f Lopt/s\n', ii, muscleNames{m}, globalMax(m));
end
fprintf('\nGroup1 (Pec+Tri):\n'); disp(muscleNames(idxGroup1));
fprintf('Group2 (Del+Others):\n'); disp(muscleNames(idxGroup2));

% theta
theta = (0:K-1) * 2*pi / K;
theta = theta(:);
theta_c = [theta; theta(1)];


tmp = maxVnorm_all(idxTop,:);
rmaxCommon = max(tmp(:));

if ~isfinite(rmaxCommon) || rmaxCommon <= 0
    rmaxCommon = 1;
end

fig1 = plotPolarGroup(theta, theta_c, maxVnorm_all, idxGroup1, muscleNames, cmapMuscles, ...
    K, rmaxCommon, 'Group 1: Pectoralis + Triceps (same color, different line styles)');

fig2 = plotPolarGroup(theta, theta_c, maxVnorm_all, idxGroup2, muscleNames, cmapMuscles, ...
    K, rmaxCommon, 'Group 2: Deltoids + Others (same color, different line styles)');

%% 7. PNG+PDF
outDir = 'fig_export';
if ~exist(outDir,'dir'); mkdir(outDir); end

exportPolarFigure(fig1, fullfile(outDir,'mouse_polar_group1_PecTri'));
exportPolarFigure(fig2, fullfile(outDir,'mouse_polar_group2_DelOthers'));


%% ===== Local Functions =====

function [Q, Qdot] = makeSCurveTrajectory(time, q_start, q_end)
N  = numel(time);
nq = numel(q_start);
Q  = zeros(N, nq);

tau = (time - time(1)) / (time(end) - time(1));
s   = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;

for k = 1:N
    Q(k,:) = (q_start + (q_end - q_start) * s(k)).';
end
Qdot = finiteDiff(time, Q);
end

function Qdot = finiteDiff(time, Q)
[N, nq] = size(Q);
Qdot = zeros(N, nq);
for j = 1:nq
    qj = Q(:,j);
    for k = 2:N-1
        dt = time(k+1) - time(k-1);
        Qdot(k,j) = (qj(k+1) - qj(k-1)) / dt;
    end
    Qdot(1,j) = (qj(2) - qj(1)) / (time(2) - time(1));
    Qdot(N,j) = (qj(N) - qj(N-1)) / (time(N) - time(N-1));
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

function fig = plotPolarGroup(theta, theta_c, maxVnorm_all, idxList, muscleNames, cmapMuscles, ...
                             K, rmaxCommon, titleStr)
fig = figure('Name',titleStr,'Units','normalized','Position',[0.08 0.10 0.68 0.78]);
pax = polaraxes(fig);
pax.Position = [0.08 0.08 0.70 0.84];
pax.FontSize = 12;
hold(pax,'on');

pax.ThetaZeroLocation = 'top';
pax.ThetaDir          = 'clockwise';

% ===== Style Line=====
pecStyles = {'-','--',':'};   % Pectoralis：
delStyles = {'-','--'};       % Deltoid：

pecColor = [];
delColor = [];
pecCnt = 0;
delCnt = 0;

h = gobjects(numel(idxList),1);

for ii = 1:numel(idxList)
    m = idxList(ii);
    name = muscleNames{m};

    r_m = maxVnorm_all(m,:).';
    r_m(~isfinite(r_m)) = 0;
    r_m(r_m < 0) = 0;
    r_c = [r_m; r_m(1)];

    % 
    col = cmapMuscles(m,:);
    ls  = '-';

    % Pectoralis：
    if startsWith(name,'Pectoralis_')
        if isempty(pecColor); pecColor = col; end
        pecCnt = pecCnt + 1;
        col = pecColor;
        ls  = pecStyles{min(pecCnt, numel(pecStyles))};

    % Deltoid：
    elseif startsWith(name,'Deltoid_')
        if isempty(delColor); delColor = col; end
        delCnt = delCnt + 1;
        col = delColor;
        ls  = delStyles{min(delCnt, numel(delStyles))};
    end

    h(ii) = polarplot(pax, theta_c, r_c, ...
        'Color', col, ...
        'LineStyle', ls, ...
        'LineWidth', 1.9);
end

rlim(pax, [0, 1.1*rmaxCommon]);

% Angle
theta_deg_all = rad2deg(theta);
step = max(1, round(K/10));      % K=100 -> 10
idxLabel = 1:step:(K-1);
pax.ThetaTick      = theta_deg_all(idxLabel);
pax.ThetaTickLabel = arrayfun(@(k) sprintf('T%d',k), idxLabel, 'UniformOutput', false);

% Radius
rt = pax.RTick;
rmaxTick = rt(end);
pax.RTick = [0, rmaxTick];
pax.RTickLabel = {'0', sprintf('%.1f L_{opt}/s', rmaxTick)};

t = title(pax, titleStr, 'Interpreter','tex');
t.FontWeight = 'bold';

legend(h, muscleNames(idxList), 'Location','eastoutside','Interpreter','none','FontSize',10);
hold(pax,'off');
end
function exportPolarFigure(h, fileBase)
% Robust export for polar plots (works if h is a figure OR an axes handle).
% Exports whole figure so the legend is included.

if isempty(h) || ~isgraphics(h)
    error('exportPolarFigure: invalid graphics handle.');
end

% Get figure handle no matter what is passed in.
if isgraphics(h, 'figure')
    fig = h;
else
    fig = ancestor(h, 'figure');
    if isempty(fig) || ~isgraphics(fig, 'figure')
        fig = gcf;
    end
end

set(fig, 'Color', 'w');

% PNG (raster)
try
    exportgraphics(fig, [fileBase '.png'], ...
        'Resolution', 600, 'BackgroundColor', 'white');
catch
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, [fileBase '.png'], '-dpng', '-r600');
end

% PDF (vector preferred)
try
    exportgraphics(fig, [fileBase '.pdf'], ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
catch
    set(fig, 'Renderer', 'painters');
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, [fileBase '.pdf'], '-dpdf', '-painters');
end
end
