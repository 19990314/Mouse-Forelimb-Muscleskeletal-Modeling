% mouse_workspace_sampling.m
clear; clc; close all;
import org.opensim.modeling.*;

%% 1.Loading Model
modelPath  = 'scaled_mouse.osim';

coordNames = { ...
    'elv_angle',        ...   % /jointset/shoulder/elv_angle
    'extension_angle',  ...   % /jointset/shoulder/extension_angle
    'rotation_angle',   ...   % /jointset/shoulder/rotation_angle
    'elbow_flex',       ...   % /jointset/humerus_ulna/elbow_flex
    'wrist_angle',      ...   % /jointset/wrist/wrist_angle
    'radius_rot'        ...   % /jointset/ulna_radius_pj/radius_rot
};

markerName = 'handm';   %  marker

fprintf('Loading model: %s\n', modelPath);
model = Model(modelPath);
state = model.initSystem();

coordSet   = model.getCoordinateSet();
markerSet  = model.getMarkerSet();
engine     = model.getSimbodyEngine();

%Check
try
    marker = markerSet.get(markerName);
catch
    error('Marker "%s" not found in model.', markerName);
end
parentFrame = marker.getParentFrame();
locInParent = marker.get_location();   % Vec3

%% Loading Range of Motion
nq    = numel(coordNames);
q_min = zeros(1,nq);
q_max = zeros(1,nq);

% Initial Posture（deg）
q_start_deg = [93.884, 17.461, 3.705, 14.099, -12.890, -9.327];
q_start     = deg2rad(q_start_deg(:));   % [6x1] rad

fprintf('\nCoordinate ranges (rad):\n');
for j = 1:nq
    coord    = coordSet.get(coordNames{j});
    q_min(j) = coord.getRangeMin();
    q_max(j) = coord.getRangeMax();
    fprintf('  %-18s  [%.3f, %.3f]\n', coordNames{j}, q_min(j), q_max(j));
end

%% 3.1 Sampling range of DOF- q_start as center
% NaN = Using Default ROM

%                 elv      ext      rot     elbow    wrist   radius
span_deg = [      50,     100,      100,     60,     NaN,      NaN  ];

fprintf('\nApplying user-defined sampling spans around q_{start}:\n');
for j = 1:nq
    if ~isnan(span_deg(j))
        half_span_rad = deg2rad(span_deg(j)) / 2;
        center        = q_start(j);

        q_user_min = center - half_span_rad;
        q_user_max = center + half_span_rad;

        % Intersected with default ROM.
        old_min = q_min(j);
        old_max = q_max(j);

        q_min(j) = max(q_min(j), q_user_min);
        q_max(j) = min(q_max(j), q_user_max);

        fprintf('  %-18s  center = %7.3f deg, span = %5.1f deg  -->  [%.3f, %.3f] rad\n', ...
            coordNames{j}, q_start_deg(j), span_deg(j), q_min(j), q_max(j));

        if q_min(j) > q_max(j)
            warning('%s: user span and ROM do not overlap! 请检查设置。', coordNames{j});
        end
    else
        fprintf('  %-18s  span = NaN -> 使用模型 ROM [%.3f, %.3f]\n', ...
            coordNames{j}, q_min(j), q_max(j));
    end
end
fprintf('\n');

%% 4.  Handm coordinates in ground  [X Y Z]
for j = 1:nq
    coord = coordSet.get(coordNames{j});
    coord.setValue(state, q_start(j));
end
model.realizePosition(state);

p0 = Vec3();
engine.getPosition(state, parentFrame, locInParent, p0);
startXYZ = [p0.get(0), p0.get(1), p0.get(2)];   % [X Y Z]

fprintf('Start posture handm in ground: [%.4f  %.4f  %.4f] m\n', ...
        startXYZ(1), startXYZ(2), startXYZ(3));

%% 5. Monte Carlo Sampling using Handm
N = 100000;   % Posture

Q_samples   = zeros(N, nq);   %  (rad)
XYZ_samples = zeros(N, 3);    % handm ground  [X Y Z]

% Random Seed
tstart = clock;
seed   = round(sum(1000*tstart));
rand('state', seed);

fprintf('Sampling %d random postures in joint space...\n', N);

for i = 1:N
    % 1) Random DOF
    for j = 1:nq
        q_ij = q_min(j) + (q_max(j) - q_min(j)) * rand();
        Q_samples(i,j) = q_ij;
        coord = coordSet.get(coordNames{j});
        coord.setValue(state, q_ij);
    end

    % 2) handm in ground
    model.realizePosition(state);
    pGround = Vec3();
    engine.getPosition(state, parentFrame, locInParent, pGround);

    X = pGround.get(0);   % forward
    Y = pGround.get(1);   % up
    Z = pGround.get(2);   % outward

    XYZ_samples(i,:) = [X, Y, Z];  
end

fprintf('Done sampling workspace.\n');

%% 6. Choose K target
K = 100;
idxTargets    = randperm(N, K);
q_end_tasks   = Q_samples(idxTargets,:);    % [K x nq] rad
xyz_end_tasks = XYZ_samples(idxTargets,:);  % [K x 3]  [X Y Z]

X = XYZ_samples(:,1);
Y = XYZ_samples(:,2);
Z = XYZ_samples(:,3);

X_end = xyz_end_tasks(:,1);
Y_end = xyz_end_tasks(:,2);
Z_end = xyz_end_tasks(:,3);

%% 7. Draw 3D workspace Cloud
figure('Name','Mouse hand workspace (OpenSim ground frame)');

scatter3(X, Y, Z, 5, Y, '.');
axis equal; grid on;
xlabel('X (forward, m)');
ylabel('Y (up, m)');
zlabel('Z (outward, m)');
title(sprintf('Workspace of marker %s (ground frame)', markerName), ...
      'Interpreter','none');
colorbar;
hold on;

%  handm 
plot3(startXYZ(1), startXYZ(2), startXYZ(3), ...
      'kp','MarkerSize',12,'MarkerFaceColor','k');

%  K target
scatter3(X_end, Y_end, Z_end, 60, 'r', 'filled');

% K target numbering
for k = 1:K
    text(X_end(k), Y_end(k), Z_end(k), sprintf(' %d', k), ...
        'Color','k', 'FontSize',14, 'HorizontalAlignment','left');
end

legend({'workspace samples','start posture','selected targets'}, 'Location','best');
view(135, 20);
hold off;

%% 8. Save
save('mouse_workspace_samples.mat', ...
     'modelPath', 'coordNames', 'markerName', ...
     'q_min', 'q_max', ...
     'q_start', 'q_start_deg', 'startXYZ', ...
     'Q_samples', 'XYZ_samples', ...
     'q_end_tasks', 'xyz_end_tasks');

fprintf('\nWorkspace samples and %d target tasks saved to mouse_workspace_samples.mat\n', K);

%% 9. Print
q_end_deg = rad2deg(q_end_tasks);   % [K x nq]

fprintf('\n===== Selected %d targets and their 6-DOF joint angles =====\n', K);
for k = 1:K
    fprintf('\nTarget %2d:\n', k);
    fprintf('  Endpoint in ground (m): [%.4f  %.4f  %.4f]\n', ...
        xyz_end_tasks(k,1), xyz_end_tasks(k,2), xyz_end_tasks(k,3));   % [X Y Z]
    fprintf('  Joint angles (deg):\n');
    for j = 1:nq
        fprintf('    %-18s = %8.3f deg\n', coordNames{j}, q_end_deg(k,j));
    end
end
