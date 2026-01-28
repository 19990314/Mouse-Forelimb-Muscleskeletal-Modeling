function plot_montecarlo_results(matfile)
% PLOT_MONTECARLO_RESULTS - Display figures from previously saved Monte Carlo results
%
% Usage:
%   plot_montecarlo_results()                           % Uses default filename
%   plot_montecarlo_results('my_results.mat')           % Specify custom file
%   plot_montecarlo_results('montecarlo_muscle_results_v2.mat')
%
% This function loads the saved results and regenerates all figures
% WITHOUT re-running the simulation.

%% ===== LOAD RESULTS =====
if nargin < 1
    matfile = 'montecarlo_muscle_results_v2.mat';
end

if ~isfile(matfile)
    error('Results file not found: %s\nRun montecarlo_muscle_work_analysis_v2() first.', matfile);
end

fprintf('Loading results from: %s\n', matfile);
load(matfile, 'results');

% Extract commonly used variables
cfg = results.config;
time = results.time;
nm = numel(cfg.muscles);
startXYZ = results.startXYZ;

% Top 5 data
nTop = 5;
top5_most_trials = results.top5_most_trials;
top5_least_trials = results.top5_least_trials;
Q_top5_most = results.Q_top5_most;
Q_top5_least = results.Q_top5_least;
xyz_top5_most = results.xyz_top5_most;
xyz_top5_least = results.xyz_top5_least;
Vm_norm_top5_most = results.Vm_norm_top5_most;
Vm_norm_top5_least = results.Vm_norm_top5_least;

% Primary representatives
most_contrast_trial = results.most_contrast_trial;
least_contrast_trial = results.least_contrast_trial;
Q_most = results.Q_most;
Q_least = results.Q_least;
Vm_norm_most = results.Vm_norm_most;
Vm_norm_least = results.Vm_norm_least;
Lm_most = results.Lm_most;

% Colors
colors = lines(nm);
red_colors = [linspace(0.4, 1, nTop)', linspace(0, 0.3, nTop)', linspace(0, 0.3, nTop)'];
blue_colors = [linspace(0, 0.3, nTop)', linspace(0, 0.3, nTop)', linspace(0.4, 1, nTop)'];

fprintf('Generating figures...\n');

%% ===== FIGURE 1: VELOCITY DYNAMICS =====
figure('Position', [100 100 1400 600], 'Name', 'Muscle Velocity Dynamics - Top 5', 'Color', 'w');

% MOST contrasting
subplot(1,2,1);
hold on;
for m = 1:nm
    plot(time, Vm_norm_most(:,m), 'Color', colors(m,:), 'LineWidth', 1.5, ...
         'DisplayName', strrep(cfg.muscles{m}, '_', ' '));
end
yline(0, 'k--', 'LineWidth', 1.5);
ylims = ylim;
patch([time(1) time(end) time(end) time(1)], [0 0 ylims(2) ylims(2)], ...
      [0.9 1 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
patch([time(1) time(end) time(end) time(1)], [ylims(1) ylims(1) 0 0], ...
      [1 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.02, 0.95, 'Eccentric', 'Units', 'normalized', 'FontSize', 12, 'Color', [0 0.6 0], 'FontWeight', 'bold');
text(0.02, 0.05, 'Concentric', 'Units', 'normalized', 'FontSize', 12, 'Color', [0.6 0 0], 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('V_{norm} (L_{opt}/s)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('MOST Contrasting (Trial %d)', most_contrast_trial), 'FontSize', 14);
legend('Location', 'eastoutside', 'FontSize', 8);
grid on;
hold off;

% LEAST contrasting
subplot(1,2,2);
hold on;
Vm_norm_least = results.Vm_norm_least;
for m = 1:nm
    plot(time, Vm_norm_least(:,m), 'Color', colors(m,:), 'LineWidth', 1.5, ...
         'DisplayName', strrep(cfg.muscles{m}, '_', ' '));
end
yline(0, 'k--', 'LineWidth', 1.5);
ylims = ylim;
patch([time(1) time(end) time(end) time(1)], [0 0 ylims(2) ylims(2)], ...
      [0.9 1 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
patch([time(1) time(end) time(end) time(1)], [ylims(1) ylims(1) 0 0], ...
      [1 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.02, 0.95, 'Eccentric', 'Units', 'normalized', 'FontSize', 12, 'Color', [0 0.6 0], 'FontWeight', 'bold');
text(0.02, 0.05, 'Concentric', 'Units', 'normalized', 'FontSize', 12, 'Color', [0.6 0 0], 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('V_{norm} (L_{opt}/s)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('LEAST Contrasting (Trial %d)', least_contrast_trial), 'FontSize', 14);
legend('Location', 'eastoutside', 'FontSize', 8);
grid on;
hold off;

%% ===== FIGURE 2: 3D JOINT CONFIGURATION =====
figure('Position', [150 150 1000 700], 'Name', '3D Joint Configuration - Top 5', 'Color', 'w');

shoulder_idx = find(contains(cfg.coordNames, 'elv_angle'), 1);
elbow_idx = find(contains(cfg.coordNames, 'elbow_flex'), 1);
wrist_idx = find(contains(cfg.coordNames, 'wrist_angle'), 1);

if ~isempty(shoulder_idx) && ~isempty(elbow_idx) && ~isempty(wrist_idx)
    hold on;
    
    % Plot TOP 5 MOST
    for i = nTop:-1:1
        Q_i = Q_top5_most{i};
        sh = rad2deg(Q_i(:, shoulder_idx));
        el = rad2deg(Q_i(:, elbow_idx));
        wr = rad2deg(Q_i(:, wrist_idx));
        lw = 1.0 + 2.0 * (nTop - i + 1) / nTop;
        plot3(sh, el, wr, '-', 'Color', red_colors(i,:), 'LineWidth', lw, ...
              'DisplayName', sprintf('MOST #%d', i));
        plot3(sh(end), el(end), wr(end), 'o', 'Color', red_colors(i,:), ...
              'MarkerSize', 8+2*(nTop-i), 'MarkerFaceColor', red_colors(i,:), 'HandleVisibility', 'off');
    end
    
    % Plot TOP 5 LEAST
    for i = nTop:-1:1
        Q_i = Q_top5_least{i};
        sh = rad2deg(Q_i(:, shoulder_idx));
        el = rad2deg(Q_i(:, elbow_idx));
        wr = rad2deg(Q_i(:, wrist_idx));
        lw = 1.0 + 2.0 * (nTop - i + 1) / nTop;
        plot3(sh, el, wr, '-', 'Color', blue_colors(i,:), 'LineWidth', lw, ...
              'DisplayName', sprintf('LEAST #%d', i));
        plot3(sh(end), el(end), wr(end), 's', 'Color', blue_colors(i,:), ...
              'MarkerSize', 8+2*(nTop-i), 'MarkerFaceColor', blue_colors(i,:), 'HandleVisibility', 'off');
    end
    
    % Mark start
    sh_start = rad2deg(Q_top5_most{1}(1, shoulder_idx));
    el_start = rad2deg(Q_top5_most{1}(1, elbow_idx));
    wr_start = rad2deg(Q_top5_most{1}(1, wrist_idx));
    plot3(sh_start, el_start, wr_start, 'kp', 'MarkerSize', 16, 'MarkerFaceColor', 'y', ...
          'LineWidth', 2, 'DisplayName', 'Start');
    
    xlabel([strrep(cfg.coordNames{shoulder_idx}, '_', ' ') ' (deg)'], 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Elbow Flex (deg)', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Wrist Angle (deg)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('3D Joint Configuration: Top %d MOST (red) vs LEAST (blue)', nTop), 'FontSize', 14);
    grid on; grid minor;
    view(135, 25);
    legend('Location', 'eastoutside', 'FontSize', 9);
    hold off;
end

%% ===== FIGURE 3: 3D HAND TRAJECTORY (HEAD-CENTERED) =====
figure('Position', [200 200 1200 800], 'Name', '3D Hand Trajectories (Head at Origin)', 'Color', 'w');

ax = subplot(1,2,1);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');

scale = 1000;  % m -> mm
head_to_hand_start = 25;  % mm

% Plot TOP 5 MOST contrasting
for i = nTop:-1:1
    xyz_rel = (xyz_top5_most{i} - startXYZ) * scale;
    xyz_head = xyz_rel + [head_to_hand_start, 0, 0];
    lw = 1.5 + 1.5 * (nTop - i + 1) / nTop;
    plot3(ax, xyz_head(:,1), xyz_head(:,2), xyz_head(:,3), '-', ...
          'Color', red_colors(i,:), 'LineWidth', lw, ...
          'DisplayName', sprintf('MOST #%d (T%d)', i, top5_most_trials(i)));
    plot3(ax, xyz_head(end,1), xyz_head(end,2), xyz_head(end,3), 'o', ...
          'Color', red_colors(i,:), 'MarkerSize', 6+2*(nTop-i), ...
          'MarkerFaceColor', red_colors(i,:), 'HandleVisibility', 'off');
end

% Plot TOP 5 LEAST contrasting
for i = nTop:-1:1
    xyz_rel = (xyz_top5_least{i} - startXYZ) * scale;
    xyz_head = xyz_rel + [head_to_hand_start, 0, 0];
    lw = 1.5 + 1.5 * (nTop - i + 1) / nTop;
    plot3(ax, xyz_head(:,1), xyz_head(:,2), xyz_head(:,3), '-', ...
          'Color', blue_colors(i,:), 'LineWidth', lw, ...
          'DisplayName', sprintf('LEAST #%d (T%d)', i, top5_least_trials(i)));
    plot3(ax, xyz_head(end,1), xyz_head(end,2), xyz_head(end,3), 's', ...
          'Color', blue_colors(i,:), 'MarkerSize', 6+2*(nTop-i), ...
          'MarkerFaceColor', blue_colors(i,:), 'HandleVisibility', 'off');
end

% HEAD at origin
plot3(ax, 0, 0, 0, 'ko', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 20, ...
      'LineWidth', 2, 'DisplayName', 'Head (origin)');
text(ax, 0, 3, 0, 'HEAD', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Hand START position
plot3(ax, head_to_hand_start, 0, 0, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 14, ...
      'LineWidth', 2, 'DisplayName', 'Hand Start');

% Body line
plot3(ax, [0, head_to_hand_start], [0, 0], [0, 0], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');

% Coordinate axes
axis_length = 10;
quiver3(ax, 0, 0, 0, axis_length, 0, 0, 0, 'Color', [0.8 0 0], 'LineWidth', 3, 'MaxHeadSize', 0.3, 'HandleVisibility', 'off');
text(ax, axis_length+2, 0, 0, '+X (Reach)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
quiver3(ax, 0, 0, 0, 0, axis_length, 0, 0, 'Color', [0 0.7 0], 'LineWidth', 3, 'MaxHeadSize', 0.3, 'HandleVisibility', 'off');
text(ax, 0, axis_length+2, 0, '+Y (Up)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0 0.7 0]);
quiver3(ax, 0, 0, 0, 0, 0, axis_length, 0, 'Color', [0 0 0.8], 'LineWidth', 3, 'MaxHeadSize', 0.3, 'HandleVisibility', 'off');
text(ax, 0, 0, axis_length+2, '+Z (Right)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0 0 0.8]);

xlabel(ax, 'X (mm) — Reaching Direction', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax, 'Y (mm) — Up/Down', 'FontSize', 12, 'FontWeight', 'bold');
zlabel(ax, 'Z (mm) — Left(-)/Right(+)', 'FontSize', 12, 'FontWeight', 'bold');
title(ax, {'Hand Trajectories (Head at Origin)'; sprintf('Top %d MOST (red) vs LEAST (blue)', nTop)}, 'FontSize', 14);
view(ax, 145, 20);
legend(ax, 'Location', 'northeast', 'FontSize', 8);
hold(ax, 'off');

% Right panel: explanation
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
};
text(ax2, 0.05, 0.88, coord_text, 'FontSize', 11, 'VerticalAlignment', 'top', ...
    'FontName', 'FixedWidth', 'Units', 'normalized');

%% ===== FIGURE 4: MUSCLE LENGTH CHANGES =====
figure('Position', [250 250 1200 500], 'Name', 'Muscle Length Changes', 'Color', 'w');

dL_most = Lm_most(end,:) - Lm_most(1,:);
[~, idxSort] = sort(abs(dL_most), 'descend');
nTopMuscles = min(5, nm);
idxTop = idxSort(1:nTopMuscles);

% ΔL(t) plot
subplot(1,2,1);
hold on;
for jj = 1:nTopMuscles
    mIdx = idxTop(jj);
    dL_t_mm = 1000 * (Lm_most(:, mIdx) - Lm_most(1, mIdx));
    plot(time, dL_t_mm, 'Color', colors(mIdx,:), 'LineWidth', 1.8, ...
         'DisplayName', cfg.muscles{mIdx});
end
yline(0, 'k--');
xlabel('Time (s)');
ylabel('\Delta L (mm)');
title(sprintf('Top %d Muscles Length Change (Trial %d)', nTopMuscles, most_contrast_trial));
legend('Location', 'eastoutside', 'Interpreter', 'none');
grid on;
hold off;

% Bar chart
subplot(1,2,2);
dL_net_mm = 1000 * dL_most;
barColors = repmat([0.7 0.7 0.7], nm, 1);
barColors(idxTop, :) = repmat([0.85 0.1 0.1], length(idxTop), 1);
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

%% ===== FIGURE 5: JOINT EXCURSION COMPARISON =====
figure('Position', [300 100 1000 500], 'Name', 'Joint Excursion Comparison', 'Color', 'w');

delta_Q_most = results.delta_Q_most;
delta_Q_least = results.delta_Q_least;

subplot(1,2,1);
bar_data = [abs(delta_Q_most); abs(delta_Q_least)]';
b = bar(bar_data);
b(1).FaceColor = [0.85 0.2 0.2];
b(2).FaceColor = [0.2 0.2 0.85];
set(gca, 'XTickLabel', strrep(cfg.coordNames, '_', ' '), 'XTickLabelRotation', 45);
ylabel('|\Delta q| (deg)', 'FontWeight', 'bold');
xlabel('Joint', 'FontWeight', 'bold');
title({'Joint Excursion Comparison'; sprintf('MOST: %.1f° total vs LEAST: %.1f° total', ...
    sum(abs(delta_Q_most)), sum(abs(delta_Q_least)))}, 'FontSize', 11);
legend('MOST', 'LEAST', 'Location', 'northeast');
grid on;

subplot(1,2,2);
axis off;
text(0.5, 0.9, 'Summary', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
summary_text = {
    ''
    sprintf('MOST Contrasting Trial: %d', most_contrast_trial)
    sprintf('  Total joint excursion: %.1f°', sum(abs(delta_Q_most)))
    ''
    sprintf('LEAST Contrasting Trial: %d', least_contrast_trial)
    sprintf('  Total joint excursion: %.1f°', sum(abs(delta_Q_least)))
    ''
    sprintf('Number of muscles analyzed: %d', nm)
    sprintf('Trajectory duration: %.2f s', time(end))
};
text(0.1, 0.75, summary_text, 'FontSize', 11, 'VerticalAlignment', 'top');

fprintf('\n✓ All figures generated from saved results.\n');
fprintf('  To save figures, use: savefig(gcf, ''filename.fig'')\n');

end