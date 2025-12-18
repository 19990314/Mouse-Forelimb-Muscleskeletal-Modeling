% mouse_param_workspace_marker.m
%
% Purpose:
%   - Use information already computed in mouse_workspace_samples.mat
%     (modelPath, coordNames, markerName, q_min, q_max)
%   - Under the same joint ROM and random joint configurations:
%       * Resample a workspace using the original model (for comparison)
%       * Create several "model parameter variants":
%           - Slightly translate the handm marker location on the hand
%             by (Δx, Δy, Δz)
%       * For each variant, randomly sample joint postures -> compute the
%         handm position in ground coordinates
%   - Plot in 3D: the original workspace + the workspaces after parameter perturbations
%
%   Changing model parameters can change the landing points / workspace.

clear; clc; close all;
import org.opensim.modeling.*;

%% 1. Load the previous workspace results to get model and DOF info
dataFile = 'mouse_workspace_samples.mat';
if ~isfile(dataFile)
    error('File "%s" not found. Please run mouse_workspace_sampling.m once first.', dataFile);
end

S = load(dataFile);

modelPath   = S.modelPath;
coordNames  = S.coordNames;
markerName  = S.markerName;     % usually 'handm'
q_min       = S.q_min;
q_max       = S.q_max;

fprintf('Loaded workspace info from %s\n', dataFile);
fprintf('  Model path : %s\n', modelPath);
fprintf('  Marker name: %s\n', markerName);
fprintf('  DOFs       : %d\n\n', numel(coordNames));

nq = numel(coordNames);

%% 2. Monte Carlo sample count (per model variant)
N_sample = 10000;   % increase/decrease depending on speed vs. accuracy

%% 3. Define a set of "model parameter variants": offsets applied to the handm marker location
%    Each row = [dx dy dz] (units: meters, expressed in the marker's parent frame on the hand)
%
%    Here are 4 examples:
%      - [ 0      0      0    ]  -> original model (no change)
%      - [ +0.002 0      0    ]  -> move handm +X by 2 mm
%      - [ 0     +0.002 0    ]  -> move handm +Y by 2 mm
%      - [ 0      0     +0.002] -> move handm +Z by 2 mm

delta_marker = [
    0       0        0;
    0.002   0        0;
    0       0.002    0;
    0       0        0.002
];

nModelVar = size(delta_marker,1);

fprintf('We will evaluate %d model variants (different marker locations).\n\n', nModelVar);

% Store the workspace point clouds for each model variant
XYZ_var = cell(nModelVar,1);

%% 4. Loop over each "model variant": load model -> modify marker parameter -> Monte Carlo sample workspace

for r = 1:nModelVar
    fprintf('=== Variant %d / %d: Δmarker = [%.3f  %.3f  %.3f] m ===\n', ...
        r, nModelVar, delta_marker(r,1), delta_marker(r,2), delta_marker(r,3));

    % 4.1 Build a fresh model from file each time (avoid contaminating other variants)
    model = Model(modelPath);
    state = model.initSystem();

    coordSet   = model.getCoordinateSet();
    markerSet  = model.getMarkerSet();
    engine     = model.getSimbodyEngine();

    % 4.2 Find the handm marker and get its original location
    try
        marker = markerSet.get(markerName);
    catch
        error('Marker "%s" not found in model.', markerName);
    end

    loc0 = marker.get_location();   % Vec3
    dx   = delta_marker(r,1);
    dy   = delta_marker(r,2);
    dz   = delta_marker(r,3);

    % 4.3 Set the new marker location (core parameter change)
    loc_new = Vec3( ...
        loc0.get(0) + dx, ...
        loc0.get(1) + dy, ...
        loc0.get(2) + dz );
    marker.set_location(loc_new);

    % Typically, after changing properties, re-initialize the system
    state = model.initSystem();

    % 4.4 Monte Carlo sample joint space, similar to mouse_workspace_sampling
    XYZ_samples = zeros(N_sample, 3);

    tstart = clock;
    seed   = round(sum(1000*tstart));
    rand('state', seed);

    fprintf('  Sampling %d random postures for this variant...\n', N_sample);

    for i = 1:N_sample
        % Generate a random q uniformly within [q_min, q_max]
        q_i = zeros(1,nq);
        for j = 1:nq
            q_ij = q_min(j) + (q_max(j) - q_min(j)) * rand();
            q_i(j) = q_ij;

            coord = coordSet.get(coordNames{j});
            coord.setValue(state, q_ij);
        end

        model.realizePosition(state);

        % handm position in ground
        parentFrame = marker.getParentFrame();
        locInParent = marker.get_location();   % now already updated

        pGround = Vec3();
        engine.getPosition(state, parentFrame, locInParent, pGround);

        X = pGround.get(0);
        Y = pGround.get(1);
        Z = pGround.get(2);

        XYZ_samples(i,:) = [X, Y, Z];
    end

    XYZ_var{r} = XYZ_samples;

    fprintf('  Done variant %d.\n\n', r);
end

%% 5. Plot 3D: original workspace (if available) + each variant workspace

figure('Name','Workspace under different model parameters (marker location)');
hold on; grid on; axis equal;

% 5.1 If mouse_workspace_samples.mat already contains original XYZ_samples, plot as background
if isfield(S, 'XYZ_samples')
    X0 = S.XYZ_samples(:,1);
    Y0 = S.XYZ_samples(:,2);
    Z0 = S.XYZ_samples(:,3);

    scatter3(X0, Y0, Z0, 3, [0.8 0.8 0.8], 'filled');  % light-gray background
    legendEntries = {'Original workspace (from file)'};
else
    legendEntries = {};
end

% 5.2 Pick colors for each model variant
colors = lines(nModelVar);

for r = 1:nModelVar
    XYZ = XYZ_var{r};
    X = XYZ(:,1); Y = XYZ(:,2); Z = XYZ(:,3);

    scatter3(X, Y, Z, 5, colors(r,:), 'filled');

    if r == 1
        descr = sprintf('Variant %d (Δ = [0 0 0])', r);
    else
        descr = sprintf('Variant %d (Δ = [%.1fmm %.1fmm %.1fmm])', ...
            r, 1000*delta_marker(r,1), 1000*delta_marker(r,2), 1000*delta_marker(r,3));
    end
    legendEntries{end+1} = descr; %#ok<SAGROW>
end

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title(sprintf('Effect of model parameter (marker "%s" location) on workspace', markerName), ...
      'Interpreter','none');

legend(legendEntries, 'Location','bestoutside');

view(135, 25);
hold off;
