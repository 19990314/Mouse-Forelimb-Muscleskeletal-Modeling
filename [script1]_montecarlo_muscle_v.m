function results = montecarlo_muscle_work_analysis_v2(selectedMuscles)
% Monte Carlo analysis identifying muscles with significant concentric/eccentric
% switching during mouse forelimb reaching trajectories
%
% IMPROVED VERSION - Uses methods from mouse_reach_10targets_one_muscle.m
%
% Usage:
%   results = montecarlo_muscle_work_analysis_v2()              % Interactive muscle selection
%   results = montecarlo_muscle_work_analysis_v2(muscleList)    % Pre-specified muscles
%
% Key improvements over original:
%   1. Uses getSimbodyEngine().getPosition() for forward kinematics (like reference)
%   2. Uses musc.getLength() and musc.getFiberVelocity() for muscle kinematics
%   3. Uses finite difference for Qdot (like reference)
%   4. Properly handles OpenSim state realization
%   5. Matches trajectory plotting style from reference
%
% COORDINATE SYSTEM (OpenSim Ground Frame for Mouse Forelimb Model):
% ===================================================================
%   This is a FORELIMB-ONLY model (no head/torso). The scapula is 
%   fixed to the OpenSim ground frame via a WeldJoint.
%
%   Based on the ground_scapula joint orientation in the .osim file,
%   the coordinate system is interpreted as:
%
%     +X : FORWARD (rostral direction, towards where the head would be)
%     +Y : UP (dorsal direction, away from the ground)
%     +Z : LATERAL/OUTWARD (away from the body midline)
%
%   IMPORTANT: The head direction (-X) is INFERRED from the model
%   orientation, NOT extracted from actual head marker data.
%   The model contains only: scapula, humerus, ulna, radius, hand, clavicle.
%
%   Hand trajectories are shown RELATIVE to the hand start position.
%
%   Joint angles:
%     - elv_angle: shoulder elevation (arm raising in sagittal plane)
%     - extension_angle: shoulder extension/flexion
%     - rotation_angle: shoulder internal/external rotation
%     - elbow_flex: elbow flexion/extension
%     - wrist_angle: wrist flexion/extension
%     - radius_rot: forearm pronation/supination
% ===================================================================

%% ===== CONFIGURATION =====
cfg.modelPath = 'scaled_mouse.osim';
cfg.N_trials = 500;      % Monte Carlo samples
cfg.T_duration = 0.3;    % seconds per reach (matches reference: T_total = 0.3)
cfg.N_time = 201;        % time points per trajectory (matches reference: N_time = 201)
cfg.markerName = 'handm'; % hand marker for forward kinematics

% Coordinate names (same order as model)
cfg.coordNames = {'elv_angle','extension_angle','rotation_angle', ...
                  'elbow_flex','wrist_angle','radius_rot'};

% Start pose (degrees) - from reference file
cfg.q_start_deg = [93.884, 17.461, 3.705, 14.099, -12.890, -9.327];

% Workspace sampling ranges: [center_deg, span_deg] for each coordinate
% These match the reference file's span_deg values
cfg.q_ranges_deg = [
    93.884,  50;    % elv_angle
    17.461, 100;    % extension_angle
     3.705, 100;    % rotation_angle
    14.099,  60;    % elbow_flex
   -12.890,  40;    % wrist_angle (was NaN in reference, using reasonable span)
    -9.327,  30     % radius_rot (was NaN in reference, using reasonable span)
];

%% ===== LOAD MODEL =====
fprintf('Loading OpenSim model: %s\n', cfg.modelPath);
import org.opensim.modeling.*
model = Model(cfg.modelPath);
state = model.initSystem();

% Get all muscles from model
allMuscles = getAllMuscles(model);
muscles = model.getMuscles();
nm_total = muscles.getSize();
fprintf('Model has %d muscles.\n', nm_total);

% Get optimal fiber lengths for ALL muscles (for normalization)
lopt_all = zeros(1, nm_total);
muscleNames_all = cell(nm_total, 1);
for m = 0:nm_total-1
    musc = muscles.get(m);
    lopt_all(m+1) = musc.getOptimalFiberLength();
    muscleNames_all{m+1} = char(musc.getName());
end

%% ===== MUSCLE SELECTION =====
if nargin == 0 || isempty(selectedMuscles)
    % Interactive selection mode
    fprintf('\n========================================\n');
    fprintf('MUSCLES IN MODEL (%d total):\n', numel(allMuscles));
    fprintf('========================================\n');
    
    % Group and display muscles by anatomical region
    [grouped, groups] = groupMusclesByRegion(allMuscles);
    
    for g = 1:numel(groups)
        fprintf('\n--- %s ---\n', groups{g});
        muscleList = grouped.(groups{g});
        for i = 1:numel(muscleList)
            fprintf('  [%2d] %s\n', i, muscleList{i});
        end
    end
    
    fprintf('\n========================================\n');
    fprintf('SELECT MUSCLES FOR ANALYSIS:\n');
    fprintf('  Option 1: Enter indices (e.g., 1,3,5-8)\n');
    fprintf('  Option 2: Enter muscle names (comma-separated)\n');
    fprintf('  Option 3: Enter group name (e.g., "Shoulder")\n');
    fprintf('  Option 4: Type "all" to analyze all muscles\n');
    fprintf('========================================\n');
    
    selection = input('Your selection: ', 's');
    cfg.muscles = parseSelection(selection, allMuscles, grouped, groups);
else
    % Use pre-specified muscles
    cfg.muscles = selectedMuscles;
end

fprintf('\n✓ Selected %d muscles for analysis:\n', numel(cfg.muscles));
for i = 1:numel(cfg.muscles)
    fprintf('  • %s\n', cfg.muscles{i});
end
fprintf('\n');

nm = numel(cfg.muscles);

% Find indices of selected muscles in the full muscle set
selectedMuscleIdx = zeros(nm, 1);
for i = 1:nm
    for j = 1:nm_total
        if strcmp(cfg.muscles{i}, muscleNames_all{j})
            selectedMuscleIdx(i) = j;
            break;
        end
    end
end

% Get lopt for selected muscles
lopt = lopt_all(selectedMuscleIdx);

%% ===== EXTRACT MUSCLE PARAMETERS =====
params = extractMuscleParams(model, cfg.muscles);

%% ===== INITIALIZE STORAGE =====
results.config = cfg;
results.params = params;
results.muscleNames = cfg.muscles;
results.lopt = lopt;
results.scores = zeros(cfg.N_trials, nm);       % versatility score
results.W_conc = zeros(cfg.N_trials, nm);       % concentric work
results.W_ecc  = zeros(cfg.N_trials, nm);       % eccentric work
results.v_range = zeros(cfg.N_trials, nm);      % velocity range
results.switch_count = zeros(cfg.N_trials, nm); % number of conc/ecc switches
results.dL_net = zeros(cfg.N_trials, nm);       % net length change
results.maxV = zeros(cfg.N_trials, nm);         % max normalized velocity
results.trajectories = cell(cfg.N_trials, 1);   % Q(t) for all trajectories
results.xyz_trajectories = cell(cfg.N_trials, 1); % hand marker XYZ

%% ===== SETUP =====
q_start = deg2rad(cfg.q_start_deg(:));  % [6 x 1] rad
time = linspace(0, cfg.T_duration, cfg.N_time)';  % [N_time x 1]

% Get starting hand position for relative trajectory plotting
coordSet = model.getCoordinateSet();
for j = 1:numel(cfg.coordNames)
    coord = coordSet.get(cfg.coordNames{j});
    coord.setValue(state, q_start(j));
end
model.realizePosition(state);

markerSet = model.getMarkerSet();
marker = markerSet.get(cfg.markerName);
engine = model.getSimbodyEngine();
parentFrame = marker.getParentFrame();
locInParent = marker.get_location();

pGround = Vec3();
engine.getPosition(state, parentFrame, locInParent, pGround);
startXYZ = [pGround.get(0), pGround.get(1), pGround.get(2)];
fprintf('Start hand position (m): [%.4f, %.4f, %.4f]\n', startXYZ(1), startXYZ(2), startXYZ(3));

results.startXYZ = startXYZ;

%% ===== MONTE CARLO LOOP =====
fprintf('\nRunning Monte Carlo simulation (%d trials)...\n', cfg.N_trials);

% Storage for quality flags
results.valid_trajectory = true(cfg.N_trials, 1);  % Track which trajectories are valid
results.spike_detected = false(cfg.N_trials, 1);   % Track spike detection

% Spike detection parameters - STRICT thresholds to ensure smooth physiological trajectories
% For a 0.3s movement with 201 time points, dt ≈ 0.0015s
% A smooth S-curve should have gradual velocity changes
MAX_VELOCITY_JUMP = 0.5;   % Maximum allowed jump in normalized velocity between consecutive time points (L_opt/s)
MAX_ACCELERATION = 30.0;   % Maximum allowed acceleration (L_opt/s^2)
SMOOTHNESS_FACTOR = 5.0;   % Jerk threshold multiplier (spike if jerk > SMOOTHNESS_FACTOR * median_jerk)

for trial = 1:cfg.N_trials
    if mod(trial, 50) == 0
        fprintf('  Trial %d/%d\n', trial, cfg.N_trials);
    end
    
    % 1. Generate random target pose within workspace
    q_end = sampleWorkspace(cfg.q_ranges_deg, q_start);
    
    % 2. Create S-curve trajectory (using reference method)
    [Q, Qdot] = makeSCurveTrajectory(time, q_start, q_end);
    
    % 3. Compute muscle kinematics for ALL muscles (like reference)
    [Lm_all, Vm_all] = computeMuscleKinematicsForTrajectory(...
        model, state, cfg.coordNames, time, Q, Qdot);
    
    % Extract selected muscles
    Lm = Lm_all(:, selectedMuscleIdx);
    Vm = Vm_all(:, selectedMuscleIdx);
    
    % Normalize velocity by optimal fiber length (like reference: Vm_norm_k = Vm_k ./ lopt)
    Vm_norm = Vm ./ lopt;  % [N_time x nm] in (L_opt/s)
    
    % ===== SPIKE DETECTION =====
    % Check for unrealistic velocity jumps in any muscle
    spike_found = false;
    dt = time(2) - time(1);
    
    for m = 1:nm
        vn = Vm_norm(:, m);
        
        % Method 1: Check for sudden velocity jumps
        velocity_jumps = abs(diff(vn));
        if any(velocity_jumps > MAX_VELOCITY_JUMP)
            spike_found = true;
            break;
        end
        
        % Method 2: Check for unrealistic acceleration
        acceleration = abs(diff(vn)) / dt;
        if any(acceleration > MAX_ACCELERATION)
            spike_found = true;
            break;
        end
        
        % Method 3: Check for non-smooth velocity (high-frequency oscillations)
        % A smooth S-curve trajectory should have smooth muscle velocities
        % Use second derivative (jerk in velocity) as indicator
        if length(vn) > 2
            velocity_jerk = abs(diff(diff(vn)));
            median_jerk = median(velocity_jerk);
            if median_jerk > 0 && any(velocity_jerk > SMOOTHNESS_FACTOR * median_jerk + 0.1)
                spike_found = true;
                break;
            end
        end
    end
    
    results.spike_detected(trial) = spike_found;
    results.valid_trajectory(trial) = ~spike_found;
    
    % 4. Compute hand trajectory (like reference)
    xyz_abs = computeEndEffTrajectoryFromMarker(...
        model, state, cfg.coordNames, Q, cfg.markerName);
    
    % Store trajectories
    results.trajectories{trial} = Q;
    results.xyz_trajectories{trial} = xyz_abs;
    
    % 5. Analyze each muscle (only if trajectory is valid)
    if ~spike_found
        for m = 1:nm
            lf = Lm(:, m);
            vf = Vm(:, m);
            vn = Vm_norm(:, m);  % normalized velocity
            
            % Compute metrics (like reference)
            dL_net = lf(end) - lf(1);
            maxV = max(abs(vn));
            
            results.dL_net(trial, m) = dL_net;
            results.maxV(trial, m) = maxV;
            results.v_range(trial, m) = range(vn);
            
            % Skip if no significant movement
            if std(lf) < 1e-10 || maxV < 0.001
                continue;
            end
            
            % Classify concentric (shortening, v<0) vs eccentric (lengthening, v>0)
            % Using more sensitive thresholds
            is_conc = vn < -0.01;  % concentric = shortening = negative velocity
            is_ecc  = vn >  0.01;  % eccentric = lengthening = positive velocity
            
            % Count phase switches
            phase = zeros(size(vn));
            phase(is_conc) = -1;
            phase(is_ecc)  =  1;
            switch_count = sum(abs(diff(phase)) > 0);
            results.switch_count(trial, m) = switch_count;
            
            % Compute work (simplified Hill model)
            a = 0.5;  % moderate activation
            F = computeMuscleForce(lf, vf, a, params(m));
            dL = [diff(lf); 0];
            W = F .* dL;  % instantaneous work
            
            W_conc_val = sum(W(is_conc), 'omitnan');
            W_ecc_val = sum(W(is_ecc), 'omitnan');
            
            results.W_conc(trial, m) = W_conc_val;
            results.W_ecc(trial, m) = W_ecc_val;
            
            % Versatility score: prioritize muscles that show BOTH ecc AND conc
            has_ecc = any(is_ecc);
            has_conc = any(is_conc);
            
            if has_ecc && has_conc
                % Muscle shows both types - this is what we want!
                time_ecc = sum(is_ecc);
                time_conc = sum(is_conc);
                phase_balance = min(time_ecc, time_conc) / (time_ecc + time_conc);
                
                % Score combines: velocity range, phase balance, switch count
                results.scores(trial, m) = range(vn) * phase_balance * (switch_count + 1) * 100;
            else
                results.scores(trial, m) = 0;
            end
        end
    end
end

% Report spike detection statistics
n_valid = sum(results.valid_trajectory);
n_spikes = sum(results.spike_detected);
fprintf('\nSpike Detection Summary:\n');
fprintf('  Valid trajectories: %d/%d (%.1f%%)\n', n_valid, cfg.N_trials, 100*n_valid/cfg.N_trials);
fprintf('  Rejected (velocity spikes): %d/%d (%.1f%%)\n', n_spikes, cfg.N_trials, 100*n_spikes/cfg.N_trials);

if n_valid < 10
    warning('Very few valid trajectories! Consider relaxing spike detection thresholds or checking model.');
end

fprintf('Monte Carlo simulation complete.\n\n');

%% ===== ANALYSIS =====
fprintf('===== RESULTS SUMMARY =====\n');

fprintf('\nDiagnostic Summary:\n');
fprintf('  Trials with any scores > 0: %d/%d\n', sum(any(results.scores > 0, 2)), cfg.N_trials);
fprintf('  Muscles with any scores > 0: %d/%d\n', sum(any(results.scores > 0, 1)), nm);

fprintf('\nPer-muscle statistics:\n');
for m = 1:nm
    n_active = sum(results.v_range(:,m) > 0.01);
    n_both = sum(results.scores(:,m) > 0);
    avg_dL = mean(abs(results.dL_net(:,m))) * 1000;  % mm
    max_vrange = max(results.v_range(:,m));
    avg_score = mean(results.scores(:,m));
    
    fprintf('  %s:\n', cfg.muscles{m});
    fprintf('    Active trials: %d/%d, Both ecc+conc: %d\n', n_active, cfg.N_trials, n_both);
    fprintf('    Avg |ΔL|: %.3f mm, Max v_range: %.3f L0/s\n', avg_dL, max_vrange);
    fprintf('    Avg score: %.3f\n', avg_score);
end

%% ===== IDENTIFY REPRESENTATIVE TRAJECTORIES (TOP 5 MOST AND LEAST) =====
trajectory_scores = sum(results.scores, 2);

% Only consider VALID trajectories (no spikes)
valid_idx = find(results.valid_trajectory);
valid_scores = trajectory_scores(valid_idx);

% Also filter for trajectories with non-zero scores (actual ecc/conc switching)
valid_nonzero_idx = valid_idx(valid_scores > 0);

nTop = 5;  % Number of top/bottom trajectories to analyze

fprintf('\n===== TRAJECTORY SELECTION =====\n');
fprintf('Total trials: %d\n', cfg.N_trials);
fprintf('Valid (no spikes): %d\n', length(valid_idx));
fprintf('Valid with ecc/conc switching: %d\n', length(valid_nonzero_idx));

if isempty(valid_nonzero_idx) || numel(valid_nonzero_idx) < nTop
    warning('Not enough valid trajectories with ecc/conc switching.');
    fprintf('Falling back to valid trajectories with largest/smallest velocity ranges.\n');
    
    if isempty(valid_idx)
        error('No valid trajectories found! All have velocity spikes. Check model or relax thresholds.');
    end
    
    v_range_total = sum(results.v_range, 2);
    valid_v_range = v_range_total(valid_idx);
    [~, sorted_local_idx] = sort(valid_v_range, 'descend');
    
    n_available = length(valid_idx);
    top5_most_trials = valid_idx(sorted_local_idx(1:min(nTop, n_available)));
    top5_least_trials = valid_idx(sorted_local_idx(max(1, end-nTop+1):end));
else
    % Sort valid non-zero trajectories by score
    valid_nonzero_scores = trajectory_scores(valid_nonzero_idx);
    [~, sorted_local_idx] = sort(valid_nonzero_scores, 'descend');
    
    top5_most_trials = valid_nonzero_idx(sorted_local_idx(1:nTop));
    top5_least_trials = valid_nonzero_idx(sorted_local_idx(max(1, end-nTop+1):end));
end

% Primary representatives (for detailed analysis)
most_contrast_trial = top5_most_trials(1);
least_contrast_trial = top5_least_trials(1);

fprintf('\n===== TOP 5 MOST CONTRASTING TRAJECTORIES (valid only) =====\n');
for i = 1:length(top5_most_trials)
    fprintf('  #%d: Trial %d (Score: %.3f)\n', i, top5_most_trials(i), trajectory_scores(top5_most_trials(i)));
end

fprintf('\n===== TOP 5 LEAST CONTRASTING TRAJECTORIES (valid only) =====\n');
for i = 1:length(top5_least_trials)
    fprintf('  #%d: Trial %d (Score: %.3f)\n', i, top5_least_trials(i), trajectory_scores(top5_least_trials(i)));
end

results.most_contrast_trial = most_contrast_trial;
results.least_contrast_trial = least_contrast_trial;
results.top5_most_trials = top5_most_trials;
results.top5_least_trials = top5_least_trials;

%% ===== RE-COMPUTE KINEMATICS FOR TOP 5 MOST AND LEAST TRAJECTORIES =====
fprintf('\nRe-computing kinematics for top 5 most and least contrasting trajectories...\n');

% Storage for top 5 trajectories
Q_top5_most = cell(nTop, 1);
Q_top5_least = cell(nTop, 1);
xyz_top5_most = cell(nTop, 1);
xyz_top5_least = cell(nTop, 1);
Vm_norm_top5_most = cell(nTop, 1);
Vm_norm_top5_least = cell(nTop, 1);
Lm_top5_most = cell(nTop, 1);
Lm_top5_least = cell(nTop, 1);

% Process TOP 5 MOST contrasting
for i = 1:nTop
    trial = top5_most_trials(i);
    Q_i = results.trajectories{trial};
    Qdot_i = finiteDiff(time, Q_i);
    
    [Lm_all_i, Vm_all_i] = computeMuscleKinematicsForTrajectory(...
        model, state, cfg.coordNames, time, Q_i, Qdot_i);
    
    Q_top5_most{i} = Q_i;
    Lm_top5_most{i} = Lm_all_i(:, selectedMuscleIdx);
    Vm_norm_top5_most{i} = Vm_all_i(:, selectedMuscleIdx) ./ lopt;
    xyz_top5_most{i} = results.xyz_trajectories{trial};
end

% Process TOP 5 LEAST contrasting
for i = 1:nTop
    trial = top5_least_trials(i);
    Q_i = results.trajectories{trial};
    Qdot_i = finiteDiff(time, Q_i);
    
    [Lm_all_i, Vm_all_i] = computeMuscleKinematicsForTrajectory(...
        model, state, cfg.coordNames, time, Q_i, Qdot_i);
    
    Q_top5_least{i} = Q_i;
    Lm_top5_least{i} = Lm_all_i(:, selectedMuscleIdx);
    Vm_norm_top5_least{i} = Vm_all_i(:, selectedMuscleIdx) ./ lopt;
    xyz_top5_least{i} = results.xyz_trajectories{trial};
end

% Keep primary representatives for backward compatibility
Q_most = Q_top5_most{1};
Q_least = Q_top5_least{1};
Qdot_most = finiteDiff(time, Q_most);
Qdot_least = finiteDiff(time, Q_least);
Vm_norm_most = Vm_norm_top5_most{1};
Vm_norm_least = Vm_norm_top5_least{1};
Lm_most = Lm_top5_most{1};
Lm_least = Lm_top5_least{1};
xyz_most = xyz_top5_most{1};
xyz_least = xyz_top5_least{1};

% Relative coordinates for plotting (mm)
xyz_most_rel = (xyz_most - startXYZ) * 1000;
xyz_least_rel = (xyz_least - startXYZ) * 1000;

%% ===== FIGURE 1: MUSCLE VELOCITY DYNAMICS =====
figure('Position', [100 100 1400 900], 'Name', 'Muscle Fascicular Velocity', 'Color', 'w');

% Color map (like reference)
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
];
if nm <= size(cmapBase,1)
    colors = cmapBase(1:nm,:);
else
    colors = lines(nm);
end

% MOST CONTRASTING - Velocity plot
subplot(2,2,1);
hold on;
for m = 1:nm
    plot(time, Vm_norm_most(:,m), 'Color', colors(m,:), 'LineWidth', 1.5);
end
yline(0, 'k--', 'LineWidth', 1.5);
ylims = ylim;
% Shade eccentric (positive) and concentric (negative) regions
patch([time(1) time(end) time(end) time(1)], [0 0 max(ylims(2),1) max(ylims(2),1)], ...
      [0.9 1 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
patch([time(1) time(end) time(end) time(1)], [min(ylims(1),-1) min(ylims(1),-1) 0 0], ...
      [1 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(time(end)*0.02, max(ylims(2),1)*0.8, 'Eccentric', 'FontSize', 10, 'Color', [0 0.6 0], 'FontWeight', 'bold');
text(time(end)*0.02, min(ylims(1),-1)*0.8, 'Concentric', 'FontSize', 10, 'Color', [0.6 0 0], 'FontWeight', 'bold');
xlabel('Time (s)');
ylabel('V_{norm} (L_{opt}/s)');
title(sprintf('MOST Contrasting (Trial %d)', most_contrast_trial));
legend(strrep(cfg.muscles, '_', ' '), 'Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 8);
grid on;
hold off;

% LEAST CONTRASTING - Velocity plot
subplot(2,2,2);
hold on;
for m = 1:nm
    plot(time, Vm_norm_least(:,m), 'Color', colors(m,:), 'LineWidth', 1.5);
end
yline(0, 'k--', 'LineWidth', 1.5);
ylims = ylim;
patch([time(1) time(end) time(end) time(1)], [0 0 max(ylims(2),1) max(ylims(2),1)], ...
      [0.9 1 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
patch([time(1) time(end) time(end) time(1)], [min(ylims(1),-1) min(ylims(1),-1) 0 0], ...
      [1 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(time(end)*0.02, max(ylims(2),1)*0.8, 'Eccentric', 'FontSize', 10, 'Color', [0 0.6 0], 'FontWeight', 'bold');
text(time(end)*0.02, min(ylims(1),-1)*0.8, 'Concentric', 'FontSize', 10, 'Color', [0.6 0 0], 'FontWeight', 'bold');
xlabel('Time (s)');
ylabel('V_{norm} (L_{opt}/s)');
title(sprintf('LEAST Contrasting (Trial %d)', least_contrast_trial));
legend(strrep(cfg.muscles, '_', ' '), 'Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 8);
grid on;
hold off;

%% ===== SAGITTAL PLANE TRAJECTORIES - TOP 5 MOST AND LEAST =====
subplot(2,2,3:4);
scale = 1000;  % m -> mm

% Color gradients for top 5
red_colors = [linspace(1, 0.6, nTop)', linspace(0.2, 0.4, nTop)', linspace(0.2, 0.4, nTop)'];  % Red gradient for MOST
blue_colors = [linspace(0.2, 0.4, nTop)', linspace(0.2, 0.4, nTop)', linspace(1, 0.6, nTop)'];  % Blue gradient for LEAST

hold on;

% Plot TOP 5 MOST contrasting (red shades, darkest = #1)
for i = nTop:-1:1  % Plot in reverse so #1 is on top
    xyz_rel = (xyz_top5_most{i} - startXYZ) * scale;
    lw = 1.0 + 1.5 * (nTop - i + 1) / nTop;  % Thicker line for higher rank
    plot(xyz_rel(:,1), xyz_rel(:,2), '-', 'Color', red_colors(i,:), 'LineWidth', lw, ...
         'DisplayName', sprintf('MOST #%d (Trial %d)', i, top5_most_trials(i)));
    % Mark endpoint
    plot(xyz_rel(end,1), xyz_rel(end,2), 'o', 'Color', red_colors(i,:), ...
         'MarkerSize', 6+2*(nTop-i), 'MarkerFaceColor', red_colors(i,:), 'HandleVisibility', 'off');
end

% Plot TOP 5 LEAST contrasting (blue shades, darkest = #1)
for i = nTop:-1:1
    xyz_rel = (xyz_top5_least{i} - startXYZ) * scale;
    lw = 1.0 + 1.5 * (nTop - i + 1) / nTop;
    plot(xyz_rel(:,1), xyz_rel(:,2), '-', 'Color', blue_colors(i,:), 'LineWidth', lw, ...
         'DisplayName', sprintf('LEAST #%d (Trial %d)', i, top5_least_trials(i)));
    % Mark endpoint
    plot(xyz_rel(end,1), xyz_rel(end,2), 's', 'Color', blue_colors(i,:), ...
         'MarkerSize', 6+2*(nTop-i), 'MarkerFaceColor', blue_colors(i,:), 'HandleVisibility', 'off');
end

% Mark common start point
plot(0, 0, 'kp', 'MarkerSize', 16, 'MarkerFaceColor', 'y', 'LineWidth', 2, 'DisplayName', 'Start (common)');

xlabel('\Delta X (mm)');
ylabel('\Delta Y (mm)');
title(sprintf('Hand Trajectories: Top %d MOST (red) vs LEAST (blue) Contrasting', nTop));
axis equal;
grid on;
legend('Location', 'eastoutside', 'FontSize', 8);
hold off;

%% ===== FIGURE 2: 3D JOINT ANGLE SPACE - TOP 5 MOST AND LEAST =====
figure('Position', [150 150 1000 700], 'Name', '3D Joint Angle Space - Top 5', 'Color', 'w');

% Find indices for shoulder, elbow, wrist
shoulder_idx = find(contains(cfg.coordNames, 'elv_angle'), 1);
elbow_idx = find(contains(cfg.coordNames, 'elbow_flex'), 1);
wrist_idx = find(contains(cfg.coordNames, 'wrist_angle'), 1);

if ~isempty(shoulder_idx) && ~isempty(elbow_idx) && ~isempty(wrist_idx)
    hold on;
    
    % Plot TOP 5 MOST contrasting (red shades)
    for i = nTop:-1:1
        Q_i = Q_top5_most{i};
        sh = rad2deg(Q_i(:, shoulder_idx));
        el = rad2deg(Q_i(:, elbow_idx));
        wr = rad2deg(Q_i(:, wrist_idx));
        
        lw = 1.0 + 2.0 * (nTop - i + 1) / nTop;
        plot3(sh, el, wr, '-', 'Color', red_colors(i,:), 'LineWidth', lw, ...
              'DisplayName', sprintf('MOST #%d', i));
        % Mark endpoint
        plot3(sh(end), el(end), wr(end), 'o', 'Color', red_colors(i,:), ...
              'MarkerSize', 8+2*(nTop-i), 'MarkerFaceColor', red_colors(i,:), 'HandleVisibility', 'off');
    end
    
    % Plot TOP 5 LEAST contrasting (blue shades)
    for i = nTop:-1:1
        Q_i = Q_top5_least{i};
        sh = rad2deg(Q_i(:, shoulder_idx));
        el = rad2deg(Q_i(:, elbow_idx));
        wr = rad2deg(Q_i(:, wrist_idx));
        
        lw = 1.0 + 2.0 * (nTop - i + 1) / nTop;
        plot3(sh, el, wr, '-', 'Color', blue_colors(i,:), 'LineWidth', lw, ...
              'DisplayName', sprintf('LEAST #%d', i));
        % Mark endpoint
        plot3(sh(end), el(end), wr(end), 's', 'Color', blue_colors(i,:), ...
              'MarkerSize', 8+2*(nTop-i), 'MarkerFaceColor', blue_colors(i,:), 'HandleVisibility', 'off');
    end
    
    % Mark common start point
    sh_start = rad2deg(Q_top5_most{1}(1, shoulder_idx));
    el_start = rad2deg(Q_top5_most{1}(1, elbow_idx));
    wr_start = rad2deg(Q_top5_most{1}(1, wrist_idx));
    plot3(sh_start, el_start, wr_start, 'kp', 'MarkerSize', 16, 'MarkerFaceColor', 'y', ...
          'LineWidth', 2, 'DisplayName', 'Start');
    
    xlabel([strrep(cfg.coordNames{shoulder_idx}, '_', ' ') ' (deg)'], 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Elbow Flex (deg)', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Wrist Angle (deg)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('3D Joint Configuration: Top %d MOST (red) vs LEAST (blue)', nTop), 'FontSize', 14);
    grid on;
    grid minor;
    view(135, 25);
    legend('Location', 'eastoutside', 'FontSize', 9);
    hold off;
else
    warning('Could not find shoulder/elbow/wrist coordinates for 3D plot');
end

%% ===== FIGURE 3: 3D HAND TRAJECTORY WITH HEAD-CENTERED COORDINATE SYSTEM =====
% Head at origin, +X = reaching direction (away from head)

figure('Position', [200 200 1200 800], 'Name', '3D Hand Trajectories (Head at Origin)', 'Color', 'w');

% Main 3D plot
ax = subplot(1,2,1);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');

scale = 1000;  % m -> mm

% In the original model frame, -X points towards the head
% To make HEAD the origin with +X = reaching direction:
% We need to FLIP the X-axis and OFFSET so head is at origin
%
% Approximate anatomy (mm):
%   head_to_shoulder ≈ 12 mm
%   shoulder_to_hand_start ≈ 15 mm (depends on pose)
%   Total: head to hand_start ≈ 27 mm in +X direction
%
% Strategy: 
%   1. Take trajectory relative to hand start (already computed)
%   2. Add offset so that hand START is at +X from head origin
%   3. Head stays at (0,0,0)

head_to_hand_start = 25;  % mm (approximate, adjustable)

% Plot TOP 5 MOST contrasting (red shades) - smooth curves
for i = nTop:-1:1
    xyz_rel = (xyz_top5_most{i} - startXYZ) * scale;  % relative to hand start
    % Shift so hand start is at (head_to_hand_start, 0, 0) from head origin
    xyz_head = xyz_rel + [head_to_hand_start, 0, 0];
    
    lw = 1.5 + 1.5 * (nTop - i + 1) / nTop;
    plot3(ax, xyz_head(:,1), xyz_head(:,2), xyz_head(:,3), '-', ...
          'Color', red_colors(i,:), 'LineWidth', lw, ...
          'DisplayName', sprintf('MOST #%d (T%d)', i, top5_most_trials(i)));
    % Mark endpoint
    plot3(ax, xyz_head(end,1), xyz_head(end,2), xyz_head(end,3), 'o', ...
          'Color', red_colors(i,:), 'MarkerSize', 6+2*(nTop-i), ...
          'MarkerFaceColor', red_colors(i,:), 'HandleVisibility', 'off');
end

% Plot TOP 5 LEAST contrasting (blue shades) - smooth curves
for i = nTop:-1:1
    xyz_rel = (xyz_top5_least{i} - startXYZ) * scale;
    xyz_head = xyz_rel + [head_to_hand_start, 0, 0];
    
    lw = 1.5 + 1.5 * (nTop - i + 1) / nTop;
    plot3(ax, xyz_head(:,1), xyz_head(:,2), xyz_head(:,3), '-', ...
          'Color', blue_colors(i,:), 'LineWidth', lw, ...
          'DisplayName', sprintf('LEAST #%d (T%d)', i, top5_least_trials(i)));
    % Mark endpoint
    plot3(ax, xyz_head(end,1), xyz_head(end,2), xyz_head(end,3), 's', ...
          'Color', blue_colors(i,:), 'MarkerSize', 6+2*(nTop-i), ...
          'MarkerFaceColor', blue_colors(i,:), 'HandleVisibility', 'off');
end

% HEAD at origin (0,0,0)
plot3(ax, 0, 0, 0, 'ko', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 20, ...
      'LineWidth', 2, 'DisplayName', 'Head (origin)');
text(ax, 0, 3, 0, 'HEAD', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Hand START position
plot3(ax, head_to_hand_start, 0, 0, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 14, ...
      'LineWidth', 2, 'DisplayName', 'Hand Start');

% Draw body/arm line from head to hand start (simplified)
plot3(ax, [0, head_to_hand_start], [0, 0], [0, 0], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');

% Coordinate axes at HEAD (origin)
axis_length = 10;  % mm

% +X axis (Reaching direction)
quiver3(ax, 0, 0, 0, axis_length, 0, 0, 0, 'Color', [0.8 0 0], 'LineWidth', 3, ...
        'MaxHeadSize', 0.3, 'HandleVisibility', 'off');
text(ax, axis_length+2, 0, 0, '+X (Reach)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.8 0 0]);

% +Y axis (Up)
quiver3(ax, 0, 0, 0, 0, axis_length, 0, 0, 'Color', [0 0.7 0], 'LineWidth', 3, ...
        'MaxHeadSize', 0.3, 'HandleVisibility', 'off');
text(ax, 0, axis_length+2, 0, '+Y (Up)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0 0.7 0]);

% +Z axis (Lateral/Right)
quiver3(ax, 0, 0, 0, 0, 0, axis_length, 0, 'Color', [0 0 0.8], 'LineWidth', 3, ...
        'MaxHeadSize', 0.3, 'HandleVisibility', 'off');
text(ax, 0, 0, axis_length+2, '+Z (Right)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0 0 0.8]);

xlabel(ax, 'X (mm) — Reaching Direction', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax, 'Y (mm) — Up/Down', 'FontSize', 12, 'FontWeight', 'bold');
zlabel(ax, 'Z (mm) — Left(-)/Right(+)', 'FontSize', 12, 'FontWeight', 'bold');
title(ax, {'Hand Trajectories (Head at Origin)'; ...
           sprintf('Top %d MOST (red) vs LEAST (blue)', nTop)}, 'FontSize', 14);
view(ax, 145, 20);
legend(ax, 'Location', 'northeast', 'FontSize', 8);
hold(ax, 'off');

% Right panel: Coordinate system explanation
ax2 = subplot(1,2,2);
axis(ax2, 'off');

text(ax2, 0.5, 0.98, 'HEAD-CENTERED COORDINATES', 'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Units', 'normalized');

coord_text = {
    ''
    'ORIGIN: Head (inferred position)'
    ''
    '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
    ''
    '  +X (Red):   REACHING DIRECTION'
    '              Away from head (rostral)'
    ''
    '  +Y (Green): UP (dorsal)'  
    '              Away from ground'
    ''
    '  +Z (Blue):  RIGHT (lateral)'
    '              Mouse''s right side'
    ''
    '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
    ''
    'LAYOUT (along +X axis):'
    ''
    '  HEAD -----> HAND START ~~~~> REACH'
    '   0 mm        ~25 mm'
    ''
    '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
    ''
    'NOTE:'
    '  • Head position is INFERRED'
    '    (model has forelimb only)'
    '  • Distance is approximate'
    '  • Trajectories show hand movement'
    '    during 0.3s reaching motion'
};

text(ax2, 0.05, 0.88, coord_text, 'FontSize', 11, 'VerticalAlignment', 'top', ...
    'FontName', 'FixedWidth', 'Units', 'normalized');

% Simple top-view schematic
axes('Position', [0.65, 0.08, 0.28, 0.30]);
hold on; axis equal; axis off;

% Draw mouse (top view)
% Head at origin
theta = linspace(0, 2*pi, 50);
head_r = 1;
fill(head_r*cos(theta), head_r*sin(theta), [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 2);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);  % head center = origin
text(0, 1.5, 'Head (0,0,0)', 'FontSize', 9, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Nose pointing +X
plot([head_r, head_r+0.5], [0, 0], 'k-', 'LineWidth', 2);

% Body extending in -X (behind head)
rectangle('Position', [-3, -0.8, 2, 1.6], 'Curvature', 0.3, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k');

% Arm reaching in +X direction
plot([0, 3], [0, 0.2], 'k-', 'LineWidth', 3);
plot(3, 0.2, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 12);
text(3.2, 0.2, 'Hand', 'FontSize', 9);

% Coordinate arrows
quiver(0, -2, 1.5, 0, 0, 'Color', [0.8 0 0], 'LineWidth', 2, 'MaxHeadSize', 0.4);
quiver(0, -2, 0, 1, 0, 'Color', [0 0 0.8], 'LineWidth', 2, 'MaxHeadSize', 0.4);
text(1.7, -2, '+X', 'Color', [0.8 0 0], 'FontWeight', 'bold');
text(0.1, -0.9, '+Z', 'Color', [0 0 0.8], 'FontWeight', 'bold');

title('Top View (mouse looking right)', 'FontSize', 10);
xlim([-4, 5]);
ylim([-3, 2.5]);
hold off;

%% ===== FIGURE 4: MUSCLE LENGTH CHANGES (like reference Figure 5.5) =====
figure('Position', [250 250 1200 500], 'Name', 'Muscle Length Changes', 'Color', 'w');

% Find top 5 muscles by |ΔL| for MOST contrasting
dL_most = Lm_most(end,:) - Lm_most(1,:);
[~, idxSort] = sort(abs(dL_most), 'descend');
nTop = min(5, nm);
idxTop = idxSort(1:nTop);

% ΔL(t) plot
subplot(1,2,1);
hold on;
for jj = 1:nTop
    mIdx = idxTop(jj);
    dL_t_mm = 1000 * (Lm_most(:, mIdx) - Lm_most(1, mIdx));
    plot(time, dL_t_mm, 'Color', colors(mIdx,:), 'LineWidth', 1.8, ...
         'DisplayName', cfg.muscles{mIdx});
end
yline(0, 'k--');
xlabel('Time (s)');
ylabel('\Delta L (mm)');
title(sprintf('Top %d Muscles Length Change (Trial %d)', nTop, most_contrast_trial));
legend('Location', 'eastoutside', 'Interpreter', 'none');
grid on;
ax = gca;
ax.YAxis.Exponent = 0;
hold off;

% Bar chart of net ΔL for all muscles
subplot(1,2,2);
dL_net_mm = 1000 * dL_most;

% Create color array - grey for all, red for top 5
barColors = repmat([0.7 0.7 0.7], nm, 1);  % grey for all
barColors(idxTop, :) = repmat([0.85 0.1 0.1], length(idxTop), 1);  % red for top 5

% Plot bars one by one with correct colors
hold on;
for m = 1:nm
    bar(m, dL_net_mm(m), 'FaceColor', barColors(m,:), 'EdgeColor', 'none', 'BarWidth', 0.8);
end
hold off;

xlabel('Muscle');
ylabel('\Delta L_{net} (mm)');
title(sprintf('Net Length Change (Trial %d)', most_contrast_trial));
set(gca, 'XTick', 1:nm, 'XTickLabel', strrep(cfg.muscles, '_', ' '), ...
    'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
grid on;
ax = gca;
ax.YAxis.Exponent = 0;

%% ===== FIGURE 5: DIAGNOSTIC - S-CURVE IN TIME DOMAIN =====
% This figure explains why the 3D joint trajectory appears linear:
% The S-curve is in the TIME domain, not the joint space!

figure('Position', [300 100 1400 800], 'Name', 'S-Curve Diagnostic: Time Domain vs Joint Space', 'Color', 'w');

% --- Row 1 Left: Joint angles vs TIME for MOST (shows the S-curve shape) ---
subplot(2,3,1);
hold on;
coord_colors = lines(6);
for j = 1:min(6, size(Q_most, 2))
    plot(time, rad2deg(Q_most(:,j)), 'Color', coord_colors(j,:), 'LineWidth', 2);
end
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Joint Angle (deg)', 'FontWeight', 'bold');
title(sprintf('MOST Contrasting (Trial %d)', most_contrast_trial), 'FontSize', 12);
legend(strrep(cfg.coordNames, '_', ' '), 'Location', 'eastoutside', 'FontSize', 8);
grid on;

% Add delta annotation
delta_most_total = sum(abs(rad2deg(Q_most(end,:) - Q_most(1,:))));
text(0.02, 0.95, sprintf('Total |\\Delta q| = %.1f°', delta_most_total), ...
    'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');
hold off;

% --- Row 2 Left: Joint angles vs TIME for LEAST ---
subplot(2,3,4);
hold on;
for j = 1:min(6, size(Q_least, 2))
    plot(time, rad2deg(Q_least(:,j)), 'Color', coord_colors(j,:), 'LineWidth', 2);
end
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Joint Angle (deg)', 'FontWeight', 'bold');
title(sprintf('LEAST Contrasting (Trial %d)', least_contrast_trial), 'FontSize', 12);
legend(strrep(cfg.coordNames, '_', ' '), 'Location', 'eastoutside', 'FontSize', 8);
grid on;

% Add delta annotation
delta_least_total = sum(abs(rad2deg(Q_least(end,:) - Q_least(1,:))));
text(0.02, 0.95, sprintf('Total |\\Delta q| = %.1f°', delta_least_total), ...
    'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');
hold off;

% --- Row 1 middle: The S-curve profile itself ---
subplot(2,3,2);
tau_plot = linspace(0, 1, 100);
s_plot = 10*tau_plot.^3 - 15*tau_plot.^4 + 6*tau_plot.^5;
s_dot = 30*tau_plot.^2 - 60*tau_plot.^3 + 30*tau_plot.^4;  % derivative

yyaxis left;
plot(tau_plot, s_plot, 'b-', 'LineWidth', 2.5);
ylabel('Position s(\tau)', 'Color', 'b', 'FontWeight', 'bold');
ylim([0 1.1]);

yyaxis right;
plot(tau_plot, s_dot, 'r--', 'LineWidth', 2);
ylabel('Velocity ds/d\tau', 'Color', 'r', 'FontWeight', 'bold');

xlabel('Normalized Time \tau = t/T', 'FontWeight', 'bold');
title({'Minimum-Jerk (S-Curve) Profile'; 's(\tau) = 10\tau^3 - 15\tau^4 + 6\tau^5'}, 'FontSize', 11);
grid on;
legend({'Position', 'Velocity'}, 'Location', 'east');

% Add annotation explaining the key insight
annotation('textbox', [0.35, 0.55, 0.15, 0.08], ...
    'String', {'All joints share', 'the SAME s(t)!'}, ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.8 0 0], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% --- Row 1 right: 3D trajectory with TIME coloring ---
subplot(2,3,3);
if ~isempty(shoulder_idx) && ~isempty(elbow_idx) && ~isempty(wrist_idx)
    sh = rad2deg(Q_most(:, shoulder_idx));
    el = rad2deg(Q_most(:, elbow_idx));
    wr = rad2deg(Q_most(:, wrist_idx));
    
    % Color by time to show speed variation along the straight path
    scatter3(sh, el, wr, 50, time, 'filled');
    colormap(gca, 'jet');
    cb = colorbar;
    cb.Label.String = 'Time (s)';
    cb.Label.FontWeight = 'bold';
    
    hold on;
    % Connect with thin line
    plot3(sh, el, wr, 'k-', 'LineWidth', 0.5);
    % Mark start and end
    plot3(sh(1), el(1), wr(1), 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'LineWidth', 2);
    plot3(sh(end), el(end), wr(end), 'rs', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    hold off;
    
    xlabel('Shoulder (deg)');
    ylabel('Elbow (deg)');
    zlabel('Wrist (deg)');
    title({'3D Joint Space (colored by TIME)'; 'Dense points = slow; Sparse = fast'}, 'FontSize', 11);
    grid on;
    view(135, 25);
end

% --- Row 2 middle: DIRECT COMPARISON - Joint excursions MOST vs LEAST ---
subplot(2,3,5);

% Calculate delta for each joint
delta_Q_most_deg = rad2deg(Q_most(end,:) - Q_most(1,:));
delta_Q_least_deg = rad2deg(Q_least(end,:) - Q_least(1,:));

% Grouped bar chart
bar_data = [abs(delta_Q_most_deg); abs(delta_Q_least_deg)]';
b = bar(bar_data, 'grouped');
b(1).FaceColor = [0.85 0.2 0.2];  % Red for MOST
b(2).FaceColor = [0.2 0.2 0.85];  % Blue for LEAST

set(gca, 'XTickLabel', strrep(cfg.coordNames, '_', ' '), 'XTickLabelRotation', 45);
ylabel('|\Delta q| (deg)', 'FontWeight', 'bold');
xlabel('Joint', 'FontWeight', 'bold');
title({'Joint Excursion Comparison'; sprintf('MOST: %.1f° total vs LEAST: %.1f° total', ...
    sum(abs(delta_Q_most_deg)), sum(abs(delta_Q_least_deg)))}, 'FontSize', 11);
legend('MOST', 'LEAST', 'Location', 'northeast');
grid on;

% --- Row 2 right: Explanation text box ---
subplot(2,3,6);
axis off;
text(0.5, 0.9, 'Why is the 3D joint trajectory a STRAIGHT LINE?', ...
    'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

explanation = {
    ''
    'The minimum-jerk trajectory is defined as:'
    ''
    '    Q(t) = Q_{start} + (Q_{end} - Q_{start}) \times s(t)'
    ''
    'where s(t) is the S-curve function (same for ALL joints).'
    ''
    'In joint space (removing time), this becomes:'
    ''
    '    Q = Q_{start} + \Delta Q \times s'
    ''
    'This is the parametric equation of a STRAIGHT LINE!'
    ''
    'The S-curve only controls HOW FAST you move along'
    'this line (slow-fast-slow), not the PATH shape.'
    ''
    '-------------------------------------------'
    'To get CURVED paths in joint space, you need:'
    '  • Via-points (intermediate waypoints)'
    '  • Different timing profiles per joint'
    '  • Task-space planning with IK'
};

text(0.05, 0.75, explanation, 'FontSize', 10, 'VerticalAlignment', 'top', ...
    'FontName', 'FixedWidth', 'Interpreter', 'tex');

%% ===== FIGURE 6: JOINT SPACE VS CARTESIAN SPACE - TOP 5 =====
% The hand trajectory in Cartesian space CAN be curved due to nonlinear FK!

figure('Position', [350 150 1400 600], 'Name', 'Joint Space vs Cartesian Space - Top 5', 'Color', 'w');

% Left: Joint space (linear paths)
subplot(1,2,1);
if ~isempty(shoulder_idx) && ~isempty(elbow_idx) && ~isempty(wrist_idx)
    hold on;
    
    % Plot TOP 5 MOST
    for i = nTop:-1:1
        Q_i = Q_top5_most{i};
        sh = rad2deg(Q_i(:, shoulder_idx));
        el = rad2deg(Q_i(:, elbow_idx));
        wr = rad2deg(Q_i(:, wrist_idx));
        lw = 1.0 + 1.5 * (nTop - i + 1) / nTop;
        plot3(sh, el, wr, '-', 'Color', red_colors(i,:), 'LineWidth', lw);
    end
    
    % Plot TOP 5 LEAST
    for i = nTop:-1:1
        Q_i = Q_top5_least{i};
        sh = rad2deg(Q_i(:, shoulder_idx));
        el = rad2deg(Q_i(:, elbow_idx));
        wr = rad2deg(Q_i(:, wrist_idx));
        lw = 1.0 + 1.5 * (nTop - i + 1) / nTop;
        plot3(sh, el, wr, '-', 'Color', blue_colors(i,:), 'LineWidth', lw);
    end
    
    % Mark start
    sh_start = rad2deg(Q_top5_most{1}(1, shoulder_idx));
    el_start = rad2deg(Q_top5_most{1}(1, elbow_idx));
    wr_start = rad2deg(Q_top5_most{1}(1, wrist_idx));
    plot3(sh_start, el_start, wr_start, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y');
    
    hold off;
    xlabel('Shoulder (deg)');
    ylabel('Elbow (deg)');
    zlabel('Wrist (deg)');
    title({'JOINT SPACE: Linear Paths'; '(All joints share same S-curve timing)'}, 'FontSize', 12);
    grid on;
    view(135, 25);
end

% Right: Cartesian space (can be curved!)
subplot(1,2,2);
hold on;

% Plot TOP 5 MOST
for i = nTop:-1:1
    xyz_rel = (xyz_top5_most{i} - startXYZ) * 1000;
    lw = 1.0 + 1.5 * (nTop - i + 1) / nTop;
    plot3(xyz_rel(:,1), xyz_rel(:,2), xyz_rel(:,3), '-', 'Color', red_colors(i,:), 'LineWidth', lw);
end

% Plot TOP 5 LEAST
for i = nTop:-1:1
    xyz_rel = (xyz_top5_least{i} - startXYZ) * 1000;
    lw = 1.0 + 1.5 * (nTop - i + 1) / nTop;
    plot3(xyz_rel(:,1), xyz_rel(:,2), xyz_rel(:,3), '-', 'Color', blue_colors(i,:), 'LineWidth', lw);
end

% Mark start
plot3(0, 0, 0, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y');

% Compute curvature metrics for primary trajectories
xyz_m = xyz_most_rel;
straight_line_pts = [linspace(xyz_m(1,1), xyz_m(end,1), size(xyz_m,1))', ...
                     linspace(xyz_m(1,2), xyz_m(end,2), size(xyz_m,1))', ...
                     linspace(xyz_m(1,3), xyz_m(end,3), size(xyz_m,1))'];
deviations = sqrt(sum((xyz_m - straight_line_pts).^2, 2));
max_deviation = max(deviations);
path_length = sum(sqrt(sum(diff(xyz_m).^2, 2)));
straight_length = norm(xyz_m(end,:) - xyz_m(1,:));
curvature_ratio = path_length / straight_length;

hold off;

xlabel('\Delta X (mm)');
ylabel('\Delta Y (mm)');
zlabel('\Delta Z (mm)');
title({'CARTESIAN SPACE: Can Be Curved!'; '(Nonlinear forward kinematics)'}, 'FontSize', 12);
grid on;
view(135, 25);
axis equal;

% Add metrics annotation
annotation('textbox', [0.55, 0.12, 0.4, 0.12], ...
    'String', {sprintf('MOST #1 - Max deviation: %.2f mm', max_deviation), ...
               sprintf('Path/Straight ratio: %.3f (1.0 = straight)', curvature_ratio)}, ...
    'FontSize', 10, 'EdgeColor', [0.5 0.5 0.5], 'BackgroundColor', 'w');

%% ===== FIGURE 7: PER-JOINT CONTRIBUTION ANALYSIS =====
figure('Position', [400 200 1000 600], 'Name', 'Per-Joint Contribution to Trajectory', 'Color', 'w');

% Compute delta for each joint
delta_Q_most = rad2deg(Q_most(end,:) - Q_most(1,:));
delta_Q_least = rad2deg(Q_least(end,:) - Q_least(1,:));

% Bar chart of joint excursions
subplot(2,2,1);
bar_data = [abs(delta_Q_most); abs(delta_Q_least)]';
b = bar(bar_data);
b(1).FaceColor = [0.85 0.2 0.2];
b(2).FaceColor = [0.2 0.2 0.85];
set(gca, 'XTickLabel', strrep(cfg.coordNames, '_', ' '), 'XTickLabelRotation', 45);
ylabel('|\Delta q| (deg)', 'FontWeight', 'bold');
title('Joint Excursions: MOST vs LEAST Contrasting', 'FontSize', 12);
legend('MOST', 'LEAST', 'Location', 'northeast');
grid on;

% Pie chart of relative contributions (MOST)
subplot(2,2,2);
contributions_most = abs(delta_Q_most) / sum(abs(delta_Q_most)) * 100;
pie(contributions_most);
legend(strrep(cfg.coordNames, '_', ' '), 'Location', 'eastoutside', 'FontSize', 8);
title({'MOST Contrasting'; 'Joint Contribution (%)'}, 'FontSize', 12);

% Pie chart of relative contributions (LEAST)
subplot(2,2,3);
contributions_least = abs(delta_Q_least) / sum(abs(delta_Q_least)) * 100;
pie(contributions_least);
legend(strrep(cfg.coordNames, '_', ' '), 'Location', 'eastoutside', 'FontSize', 8);
title({'LEAST Contrasting'; 'Joint Contribution (%)'}, 'FontSize', 12);

% Text summary
subplot(2,2,4);
axis off;

text(0.5, 0.95, 'Joint Contribution Analysis', ...
    'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Find dominant joints
[~, sort_idx_most] = sort(abs(delta_Q_most), 'descend');
[~, sort_idx_least] = sort(abs(delta_Q_least), 'descend');

summary_text = {
    ''
    'MOST Contrasting Trajectory:'
    sprintf('  Total excursion: %.1f deg', sum(abs(delta_Q_most)))
    sprintf('  Dominant joint: %s (%.1f%%)', cfg.coordNames{sort_idx_most(1)}, contributions_most(sort_idx_most(1)))
    sprintf('  2nd joint: %s (%.1f%%)', cfg.coordNames{sort_idx_most(2)}, contributions_most(sort_idx_most(2)))
    ''
    'LEAST Contrasting Trajectory:'
    sprintf('  Total excursion: %.1f deg', sum(abs(delta_Q_least)))
    sprintf('  Dominant joint: %s (%.1f%%)', cfg.coordNames{sort_idx_least(1)}, contributions_least(sort_idx_least(1)))
    sprintf('  2nd joint: %s (%.1f%%)', cfg.coordNames{sort_idx_least(2)}, contributions_least(sort_idx_least(2)))
    ''
    '-------------------------------------------'
    'If one joint dominates (>70%), the 3D plot'
    'will appear nearly 1D (along one axis).'
    ''
    'If two joints dominate, the 3D plot will'
    'appear nearly 2D (in a plane).'
};

text(0.05, 0.8, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', 'FontName', 'FixedWidth');

% Check for dominance
if contributions_most(sort_idx_most(1)) > 70
    text(0.5, 0.08, sprintf('⚠ WARNING: %s dominates MOST trajectory!', ...
        strrep(cfg.coordNames{sort_idx_most(1)}, '_', ' ')), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.8 0.4 0], 'HorizontalAlignment', 'center');
end

%% ===== SAVE RESULTS =====
results.time = time;
results.Q_most = Q_most;
results.Q_least = Q_least;
results.Qdot_most = Qdot_most;
results.Qdot_least = Qdot_least;
results.xyz_most = xyz_most;
results.xyz_least = xyz_least;
results.Vm_norm_most = Vm_norm_most;
results.Vm_norm_least = Vm_norm_least;
results.Lm_most = Lm_most;
results.Lm_least = Lm_least;
results.delta_Q_most = delta_Q_most;
results.delta_Q_least = delta_Q_least;

% Top 5 trajectory data
results.Q_top5_most = Q_top5_most;
results.Q_top5_least = Q_top5_least;
results.xyz_top5_most = xyz_top5_most;
results.xyz_top5_least = xyz_top5_least;
results.Vm_norm_top5_most = Vm_norm_top5_most;
results.Vm_norm_top5_least = Vm_norm_top5_least;

save('montecarlo_muscle_results_v2.mat', 'results');
fprintf('\n✓ Results saved to montecarlo_muscle_results_v2.mat\n');

end

%% ========================================================================
% LOCAL FUNCTIONS (matching reference file style)
% ========================================================================

function q_end = sampleWorkspace(ranges_deg, q_start)
    % Sample random endpoint within specified ranges
    % ranges_deg: [center, span] for each coordinate
    nq = size(ranges_deg, 1);
    q_end = zeros(nq, 1);
    for i = 1:nq
        center = deg2rad(ranges_deg(i,1));
        span = deg2rad(ranges_deg(i,2));
        q_end(i) = center + span * (rand() - 0.5);
    end
end

function [Q, Qdot] = makeSCurveTrajectory(time, q_start, q_end)
    % Minimum-jerk (S-curve) trajectory - MATCHES REFERENCE EXACTLY
    N  = numel(time);
    nq = numel(q_start);
    Q  = zeros(N, nq);
    
    tau = (time - time(1)) / (time(end) - time(1));   % 0->1
    s   = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;           % quintic polynomial
    
    for k = 1:N
        Q(k,:) = (q_start + (q_end - q_start) * s(k)).';
    end
    
    Qdot = finiteDiff(time, Q);
end

function Qdot = finiteDiff(time, Q)
    % Central finite difference - MATCHES REFERENCE EXACTLY
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
    % Compute muscle kinematics - MATCHES REFERENCE EXACTLY
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
    % Compute end-effector trajectory - MATCHES REFERENCE EXACTLY
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

function params = extractMuscleParams(model, muscleNames)
    % Extract muscle parameters from OpenSim model
    muscles = model.getMuscles();
    nm = numel(muscleNames);
    params = struct([]);
    
    for i = 1:nm
        try
            m = muscles.get(muscleNames{i});
            params(i).name = muscleNames{i};
            params(i).Fmax = m.getMaxIsometricForce();
            params(i).l0 = m.getOptimalFiberLength();
            params(i).vmax = m.getMaxContractionVelocity();
            params(i).phi0 = m.getPennationAngleAtOptimalFiberLength();
        catch
            warning('Muscle %s not found, using defaults', muscleNames{i});
            params(i).name = muscleNames{i};
            params(i).Fmax = 1.0;
            params(i).l0 = 0.01;
            params(i).vmax = 10;
            params(i).phi0 = 0;
        end
    end
end

function F = computeMuscleForce(lf, vf, a, param)
    % Simplified Hill-type force model
    ln = lf ./ param.l0;
    vn = vf ./ (param.vmax * param.l0);
    
    % Force-length (Gaussian)
    w = 0.25;
    fl = exp(-((ln - 1)./w).^2);
    
    % Force-velocity (simplified Hill)
    fv = max(0.1, min(1.8, (1 - 0.3*vn)./(1 + 1.5*vn)));
    
    % Passive (exponential beyond optimal)
    fp = (ln > 1) .* (exp(5*(ln-1)) - 1) * 0.03;
    
    F = param.Fmax * (a .* fl .* fv + fp);
end

%% ===== MUSCLE LISTING HELPERS =====

function allMuscles = getAllMuscles(model)
    muscles = model.getMuscles();
    nm = muscles.getSize();
    allMuscles = cell(nm, 1);
    for i = 0:nm-1
        allMuscles{i+1} = char(muscles.get(i).getName());
    end
end

function [grouped, groupNames] = groupMusclesByRegion(muscleNames)
    grouped = struct();
    patterns = {
        'Shoulder', {'Deltoid', 'Infraspinatus', 'Subscapularis', 'Supraspinatus'};
        'Elbow_Flexors', {'Biceps', 'Brachialis', 'Brachioradialis'};
        'Elbow_Extensors', {'Triceps', 'Anconeus'};
        'Pectoral', {'Pectoralis', 'Pec'};
        'Back', {'Latissimus', 'Trapezius', 'Rhomboid'};
        'Wrist_Forearm', {'Flexor', 'Extensor', 'Pronator', 'Supinator', 'Carpi', 'Digitorum'};
        'Other', {};
    };
    
    for p = 1:size(patterns, 1)
        grouped.(patterns{p,1}) = {};
    end
    
    for i = 1:numel(muscleNames)
        name = muscleNames{i};
        assigned = false;
        for p = 1:size(patterns, 1)-1
            keywords = patterns{p,2};
            for k = 1:numel(keywords)
                if contains(name, keywords{k}, 'IgnoreCase', true)
                    grouped.(patterns{p,1}){end+1} = name;
                    assigned = true;
                    break;
                end
            end
            if assigned, break; end
        end
        if ~assigned
            grouped.Other{end+1} = name;
        end
    end
    
    groupNames = fieldnames(grouped);
    toRemove = false(size(groupNames));
    for i = 1:numel(groupNames)
        if isempty(grouped.(groupNames{i}))
            toRemove(i) = true;
        end
    end
    groupNames(toRemove) = [];
end

function selected = parseSelection(input, allMuscles, grouped, groupNames)
    input = strtrim(input);
    
    if strcmpi(input, 'all')
        selected = allMuscles;
        return;
    end
    
    for g = 1:numel(groupNames)
        if strcmpi(input, groupNames{g})
            selected = grouped.(groupNames{g});
            return;
        end
    end
    
    if contains(input, ',') && ~contains(input, '-') && isempty(str2num(input))
        parts = strsplit(input, ',');
        selected = {};
        for i = 1:numel(parts)
            name = strtrim(parts{i});
            matches = find(contains(allMuscles, name, 'IgnoreCase', true));
            if ~isempty(matches)
                selected{end+1} = allMuscles{matches(1)};
            else
                warning('Muscle "%s" not found, skipping', name);
            end
        end
        return;
    end
    
    try
        indices = parseIndices(input, numel(allMuscles));
        selected = allMuscles(indices);
    catch
        error('Invalid selection format. Please try again.');
    end
end

function indices = parseIndices(str, maxVal)
    parts = strsplit(str, ',');
    indices = [];
    for i = 1:numel(parts)
        part = strtrim(parts{i});
        if contains(part, '-')
            range_parts = strsplit(part, '-');
            start = str2double(range_parts{1});
            stop = str2double(range_parts{2});
            indices = [indices, start:stop];
        else
            indices = [indices, str2double(part)];
        end
    end
    indices = unique(indices);
    indices(indices < 1 | indices > maxVal) = [];
end