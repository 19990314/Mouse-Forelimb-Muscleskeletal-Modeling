function export_trajectories_to_stl(matfile, options)
% EXPORT_TRAJECTORIES_TO_STL - Convert hand trajectories to 3D printable STL files
%
% Usage:
%   export_trajectories_to_stl()                          % Use default settings
%   export_trajectories_to_stl('my_results.mat')          % Specify results file
%   export_trajectories_to_stl('results.mat', options)    % Custom options
%
% Options (struct):
%   .tube_radius    - Radius of trajectory tubes (mm), default: 0.5
%   .n_sides        - Number of sides for tube cross-section, default: 12
%   .scale_factor   - Scale factor for the model, default: 1.0
%   .output_dir     - Output directory for STL files, default: 'stl_export'
%   .combine_all    - Export combined STL with all trajectories, default: true
%   .head_radius    - Radius of head sphere (mm), default: 3.0
%   .add_base       - Add a base plate for stability, default: true
%
% Output:
%   - Individual STL files for each trajectory (MOST_1.stl, LEAST_1.stl, etc.)
%   - Combined STL file with all trajectories (combined_trajectories.stl)
%   - Head marker STL (head_marker.stl)
%   - Optional base plate (base_plate.stl)
%
% After export, import STL files into OnShape, Cura, or other slicer software
% to combine and prepare for 3D printing.
%
% Required: surf2stl (download from MATLAB File Exchange if not available)
%           https://www.mathworks.com/matlabcentral/fileexchange/4512-surf2stl

%% ===== PARSE INPUTS =====
if nargin < 1 || isempty(matfile)
    matfile = 'montecarlo_muscle_results_v2.mat';
end

if nargin < 2
    options = struct();
end

% Default options
defaults.tube_radius = 0.5;       % mm
defaults.n_sides = 12;            % polygon sides for tube
defaults.scale_factor = 1.0;      % scale multiplier
defaults.output_dir = 'stl_export';
defaults.combine_all = true;
defaults.head_radius = 3.0;       % mm
defaults.add_base = true;
defaults.base_thickness = 2.0;    % mm
defaults.base_margin = 5.0;       % mm beyond trajectory extent

% Merge with user options
opt = defaults;
if isstruct(options)
    fnames = fieldnames(options);
    for i = 1:length(fnames)
        opt.(fnames{i}) = options.(fnames{i});
    end
end

%% ===== LOAD RESULTS =====
if ~isfile(matfile)
    error('Results file not found: %s\nRun montecarlo_muscle_work_analysis_v2() first.', matfile);
end

fprintf('Loading results from: %s\n', matfile);
load(matfile, 'results');

startXYZ = results.startXYZ;
xyz_top5_most = results.xyz_top5_most;
xyz_top5_least = results.xyz_top5_least;
top5_most_trials = results.top5_most_trials;
top5_least_trials = results.top5_least_trials;

nTop = length(xyz_top5_most);

%% ===== CREATE OUTPUT DIRECTORY =====
if ~exist(opt.output_dir, 'dir')
    mkdir(opt.output_dir);
end
fprintf('Output directory: %s\n', opt.output_dir);

%% ===== CONVERT TRAJECTORIES TO HEAD-CENTERED COORDINATES =====
scale = 1000;  % m -> mm
head_to_hand_start = 25;  % mm (matches plotting code)

% Store all vertices for combined export
all_vertices = {};
all_faces = {};
vertex_offset = 0;

fprintf('\nProcessing trajectories...\n');

%% ===== EXPORT INDIVIDUAL TRAJECTORY STLs =====

% Process TOP 5 MOST contrasting (red)
for i = 1:nTop
    xyz_rel = (xyz_top5_most{i} - startXYZ) * scale * opt.scale_factor;
    xyz_head = xyz_rel + [head_to_hand_start * opt.scale_factor, 0, 0];
    
    % Generate tube mesh
    [vertices, faces] = trajectory_to_tube(xyz_head, opt.tube_radius * opt.scale_factor, opt.n_sides);
    
    % Export individual STL
    filename = fullfile(opt.output_dir, sprintf('MOST_%d_trial%d.stl', i, top5_most_trials(i)));
    write_stl_binary(filename, vertices, faces, sprintf('MOST_trajectory_%d', i));
    fprintf('  Exported: %s\n', filename);
    
    % Store for combined export
    if opt.combine_all
        all_vertices{end+1} = vertices;
        all_faces{end+1} = faces + vertex_offset;
        vertex_offset = vertex_offset + size(vertices, 1);
    end
end

% Process TOP 5 LEAST contrasting (blue)
for i = 1:nTop
    xyz_rel = (xyz_top5_least{i} - startXYZ) * scale * opt.scale_factor;
    xyz_head = xyz_rel + [head_to_hand_start * opt.scale_factor, 0, 0];
    
    % Generate tube mesh
    [vertices, faces] = trajectory_to_tube(xyz_head, opt.tube_radius * opt.scale_factor * 0.8, opt.n_sides);
    
    % Export individual STL
    filename = fullfile(opt.output_dir, sprintf('LEAST_%d_trial%d.stl', i, top5_least_trials(i)));
    write_stl_binary(filename, vertices, faces, sprintf('LEAST_trajectory_%d', i));
    fprintf('  Exported: %s\n', filename);
    
    % Store for combined export
    if opt.combine_all
        all_vertices{end+1} = vertices;
        all_faces{end+1} = faces + vertex_offset;
        vertex_offset = vertex_offset + size(vertices, 1);
    end
end

%% ===== CREATE HEAD MARKER (SPHERE) =====
[head_v, head_f] = create_sphere([0, 0, 0], opt.head_radius * opt.scale_factor, 20);
filename = fullfile(opt.output_dir, 'head_marker.stl');
write_stl_binary(filename, head_v, head_f, 'head_marker');
fprintf('  Exported: %s\n', filename);

if opt.combine_all
    all_vertices{end+1} = head_v;
    all_faces{end+1} = head_f + vertex_offset;
    vertex_offset = vertex_offset + size(head_v, 1);
end

%% ===== CREATE HAND START MARKER (SMALL SPHERE) =====
hand_start_pos = [head_to_hand_start * opt.scale_factor, 0, 0];
[hand_v, hand_f] = create_sphere(hand_start_pos, opt.tube_radius * opt.scale_factor * 2, 16);
filename = fullfile(opt.output_dir, 'hand_start_marker.stl');
write_stl_binary(filename, hand_v, hand_f, 'hand_start');
fprintf('  Exported: %s\n', filename);

if opt.combine_all
    all_vertices{end+1} = hand_v;
    all_faces{end+1} = hand_f + vertex_offset;
    vertex_offset = vertex_offset + size(hand_v, 1);
end

%% ===== CREATE BODY LINE (CYLINDER FROM HEAD TO HAND START) =====
[body_v, body_f] = create_cylinder([0,0,0], hand_start_pos, opt.tube_radius * opt.scale_factor * 1.5, opt.n_sides);
filename = fullfile(opt.output_dir, 'body_line.stl');
write_stl_binary(filename, body_v, body_f, 'body_line');
fprintf('  Exported: %s\n', filename);

if opt.combine_all
    all_vertices{end+1} = body_v;
    all_faces{end+1} = body_f + vertex_offset;
    vertex_offset = vertex_offset + size(body_v, 1);
end

%% ===== CREATE COORDINATE AXES =====
axis_length = 10 * opt.scale_factor;
axis_radius = opt.tube_radius * opt.scale_factor * 0.6;

% X-axis (reaching direction)
[ax_v, ax_f] = create_cylinder([0,0,0], [axis_length,0,0], axis_radius, 8);
filename = fullfile(opt.output_dir, 'axis_X.stl');
write_stl_binary(filename, ax_v, ax_f, 'axis_X');
fprintf('  Exported: %s\n', filename);
if opt.combine_all
    all_vertices{end+1} = ax_v;
    all_faces{end+1} = ax_f + vertex_offset;
    vertex_offset = vertex_offset + size(ax_v, 1);
end

% Y-axis (up)
[ay_v, ay_f] = create_cylinder([0,0,0], [0,axis_length,0], axis_radius, 8);
filename = fullfile(opt.output_dir, 'axis_Y.stl');
write_stl_binary(filename, ay_v, ay_f, 'axis_Y');
fprintf('  Exported: %s\n', filename);
if opt.combine_all
    all_vertices{end+1} = ay_v;
    all_faces{end+1} = ay_f + vertex_offset;
    vertex_offset = vertex_offset + size(ay_v, 1);
end

% Z-axis (lateral)
[az_v, az_f] = create_cylinder([0,0,0], [0,0,axis_length], axis_radius, 8);
filename = fullfile(opt.output_dir, 'axis_Z.stl');
write_stl_binary(filename, az_v, az_f, 'axis_Z');
fprintf('  Exported: %s\n', filename);
if opt.combine_all
    all_vertices{end+1} = az_v;
    all_faces{end+1} = az_f + vertex_offset;
    vertex_offset = vertex_offset + size(az_v, 1);
end

%% ===== CREATE BASE PLATE (OPTIONAL) =====
if opt.add_base
    % Find extent of all trajectories
    all_xyz = [];
    for i = 1:nTop
        xyz_rel = (xyz_top5_most{i} - startXYZ) * scale * opt.scale_factor;
        xyz_head = xyz_rel + [head_to_hand_start * opt.scale_factor, 0, 0];
        all_xyz = [all_xyz; xyz_head];
        
        xyz_rel = (xyz_top5_least{i} - startXYZ) * scale * opt.scale_factor;
        xyz_head = xyz_rel + [head_to_hand_start * opt.scale_factor, 0, 0];
        all_xyz = [all_xyz; xyz_head];
    end
    
    x_min = min(all_xyz(:,1)) - opt.base_margin;
    x_max = max(all_xyz(:,1)) + opt.base_margin;
    z_min = min(all_xyz(:,3)) - opt.base_margin;
    z_max = max(all_xyz(:,3)) + opt.base_margin;
    y_base = min(all_xyz(:,2)) - opt.base_margin;  % Below lowest point
    
    [base_v, base_f] = create_box([x_min, y_base - opt.base_thickness, z_min], ...
                                   [x_max, y_base, z_max]);
    filename = fullfile(opt.output_dir, 'base_plate.stl');
    write_stl_binary(filename, base_v, base_f, 'base_plate');
    fprintf('  Exported: %s\n', filename);
    
    if opt.combine_all
        all_vertices{end+1} = base_v;
        all_faces{end+1} = base_f + vertex_offset;
        vertex_offset = vertex_offset + size(base_v, 1);
    end
end

%% ===== EXPORT COMBINED STL =====
if opt.combine_all
    fprintf('\nCombining all meshes...\n');
    
    % Concatenate all vertices and faces
    combined_vertices = vertcat(all_vertices{:});
    combined_faces = vertcat(all_faces{:});
    
    filename = fullfile(opt.output_dir, 'combined_trajectories.stl');
    write_stl_binary(filename, combined_vertices, combined_faces, 'combined_model');
    fprintf('  Exported: %s\n', filename);
end

%% ===== SUMMARY =====
fprintf('\n========================================\n');
fprintf('STL EXPORT COMPLETE\n');
fprintf('========================================\n');
fprintf('Output directory: %s\n', opt.output_dir);
fprintf('Files exported:\n');
files = dir(fullfile(opt.output_dir, '*.stl'));
for i = 1:length(files)
    fprintf('  • %s (%.1f KB)\n', files(i).name, files(i).bytes/1024);
end

fprintf('\n========================================\n');
fprintf('NEXT STEPS FOR 3D PRINTING:\n');
fprintf('========================================\n');
fprintf('Option 1: OnShape (free online CAD)\n');
fprintf('  1. Go to onshape.com and create account\n');
fprintf('  2. Create new document\n');
fprintf('  3. Right-click in Parts list → Import\n');
fprintf('  4. Upload combined_trajectories.stl\n');
fprintf('  5. Or import individual STLs and arrange\n');
fprintf('  6. Export as STL for printing\n');
fprintf('\n');
fprintf('Option 2: Cura (free slicer)\n');
fprintf('  1. Open Ultimaker Cura\n');
fprintf('  2. File → Open File(s)\n');
fprintf('  3. Select combined_trajectories.stl\n');
fprintf('  4. Adjust scale if needed (model is in mm)\n');
fprintf('  5. Slice and print!\n');
fprintf('\n');
fprintf('Option 3: Import individual STLs\n');
fprintf('  - MOST_*.stl: High muscle switching trajectories (print in RED)\n');
fprintf('  - LEAST_*.stl: Low muscle switching trajectories (print in BLUE)\n');
fprintf('  - head_marker.stl: Origin/head position\n');
fprintf('  - axis_*.stl: Coordinate axes\n');
fprintf('  - base_plate.stl: Stability base\n');
fprintf('\n');
fprintf('Scale: 1 unit = 1 mm (model is already in mm)\n');
fprintf('Suggested print scale: 2x-5x for better visibility\n');
fprintf('========================================\n');

%% ===== OPTIONAL: PREVIEW PLOT =====
figure('Position', [100 100 800 600], 'Name', 'STL Export Preview', 'Color', 'w');
hold on;

% Plot trajectories
for i = 1:nTop
    xyz_rel = (xyz_top5_most{i} - startXYZ) * scale * opt.scale_factor;
    xyz_head = xyz_rel + [head_to_hand_start * opt.scale_factor, 0, 0];
    plot3(xyz_head(:,1), xyz_head(:,2), xyz_head(:,3), 'r-', 'LineWidth', 2);
    
    xyz_rel = (xyz_top5_least{i} - startXYZ) * scale * opt.scale_factor;
    xyz_head = xyz_rel + [head_to_hand_start * opt.scale_factor, 0, 0];
    plot3(xyz_head(:,1), xyz_head(:,2), xyz_head(:,3), 'b-', 'LineWidth', 2);
end

% Plot markers
plot3(0, 0, 0, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', [0.3 0.3 0.3]);
plot3(hand_start_pos(1), hand_start_pos(2), hand_start_pos(3), 'kp', ...
      'MarkerSize', 12, 'MarkerFaceColor', 'y');

% Plot axes
quiver3(0,0,0, axis_length,0,0, 0, 'r', 'LineWidth', 2);
quiver3(0,0,0, 0,axis_length,0, 0, 'g', 'LineWidth', 2);
quiver3(0,0,0, 0,0,axis_length, 0, 'b', 'LineWidth', 2);

xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title('STL Export Preview (Red=MOST, Blue=LEAST)');
axis equal; grid on;
view(145, 20);
legend('MOST trajectories', 'LEAST trajectories', 'Head', 'Hand Start', 'Location', 'best');
hold off;

end

%% ========================================================================
% HELPER FUNCTIONS
% ========================================================================

function [vertices, faces] = trajectory_to_tube(xyz, radius, n_sides)
% Convert a 3D trajectory to a tube mesh
% xyz: Nx3 array of points along trajectory
% radius: tube radius
% n_sides: number of sides for tube cross-section

n_points = size(xyz, 1);

% Smooth the trajectory slightly to avoid sharp corners
if n_points > 5
    xyz = smoothdata(xyz, 'gaussian', 5);
end

% Generate tube vertices
theta = linspace(0, 2*pi, n_sides+1);
theta = theta(1:end-1);  % Remove duplicate

vertices = zeros(n_points * n_sides, 3);
faces = [];

for i = 1:n_points
    % Calculate tangent direction
    if i == 1
        tangent = xyz(2,:) - xyz(1,:);
    elseif i == n_points
        tangent = xyz(end,:) - xyz(end-1,:);
    else
        tangent = xyz(i+1,:) - xyz(i-1,:);
    end
    tangent = tangent / (norm(tangent) + eps);
    
    % Find perpendicular vectors (Frenet frame approximation)
    if abs(tangent(1)) < 0.9
        ref = [1, 0, 0];
    else
        ref = [0, 1, 0];
    end
    normal = cross(tangent, ref);
    normal = normal / (norm(normal) + eps);
    binormal = cross(tangent, normal);
    
    % Generate circle of vertices at this point
    for j = 1:n_sides
        offset = radius * (cos(theta(j)) * normal + sin(theta(j)) * binormal);
        idx = (i-1) * n_sides + j;
        vertices(idx, :) = xyz(i,:) + offset;
    end
    
    % Generate faces connecting to previous ring
    if i > 1
        for j = 1:n_sides
            j_next = mod(j, n_sides) + 1;
            
            idx_curr = (i-1) * n_sides + j;
            idx_curr_next = (i-1) * n_sides + j_next;
            idx_prev = (i-2) * n_sides + j;
            idx_prev_next = (i-2) * n_sides + j_next;
            
            % Two triangles per quad
            faces = [faces; idx_prev, idx_curr, idx_prev_next];
            faces = [faces; idx_prev_next, idx_curr, idx_curr_next];
        end
    end
end

% Cap the ends
% Start cap
center_start = mean(vertices(1:n_sides, :), 1);
vertices = [vertices; center_start];
center_idx_start = size(vertices, 1);
for j = 1:n_sides
    j_next = mod(j, n_sides) + 1;
    faces = [faces; center_idx_start, j_next, j];
end

% End cap
center_end = mean(vertices((n_points-1)*n_sides+1 : n_points*n_sides, :), 1);
vertices = [vertices; center_end];
center_idx_end = size(vertices, 1);
for j = 1:n_sides
    j_next = mod(j, n_sides) + 1;
    idx_curr = (n_points-1) * n_sides + j;
    idx_next = (n_points-1) * n_sides + j_next;
    faces = [faces; center_idx_end, idx_curr, idx_next];
end

end

function [vertices, faces] = create_sphere(center, radius, n_segments)
% Create a sphere mesh
[X, Y, Z] = sphere(n_segments);
X = X * radius + center(1);
Y = Y * radius + center(2);
Z = Z * radius + center(3);

% Convert surf to vertices and faces
[faces, vertices] = surf2patch(X, Y, Z, 'triangles');
end

function [vertices, faces] = create_cylinder(p1, p2, radius, n_sides)
% Create a cylinder mesh between two points

direction = p2 - p1;
length = norm(direction);
if length < eps
    vertices = [];
    faces = [];
    return;
end
direction = direction / length;

% Find perpendicular vectors
if abs(direction(1)) < 0.9
    ref = [1, 0, 0];
else
    ref = [0, 1, 0];
end
perp1 = cross(direction, ref);
perp1 = perp1 / norm(perp1);
perp2 = cross(direction, perp1);

% Generate vertices
theta = linspace(0, 2*pi, n_sides+1);
theta = theta(1:end-1);

vertices = zeros(2 * n_sides + 2, 3);

% Bottom ring
for j = 1:n_sides
    offset = radius * (cos(theta(j)) * perp1 + sin(theta(j)) * perp2);
    vertices(j, :) = p1 + offset;
end

% Top ring
for j = 1:n_sides
    offset = radius * (cos(theta(j)) * perp1 + sin(theta(j)) * perp2);
    vertices(n_sides + j, :) = p2 + offset;
end

% Center points for caps
vertices(2*n_sides + 1, :) = p1;
vertices(2*n_sides + 2, :) = p2;

% Generate faces
faces = [];

% Side faces
for j = 1:n_sides
    j_next = mod(j, n_sides) + 1;
    
    idx_bot = j;
    idx_bot_next = j_next;
    idx_top = n_sides + j;
    idx_top_next = n_sides + j_next;
    
    faces = [faces; idx_bot, idx_top, idx_bot_next];
    faces = [faces; idx_bot_next, idx_top, idx_top_next];
end

% Bottom cap
center_bot = 2*n_sides + 1;
for j = 1:n_sides
    j_next = mod(j, n_sides) + 1;
    faces = [faces; center_bot, j_next, j];
end

% Top cap
center_top = 2*n_sides + 2;
for j = 1:n_sides
    j_next = mod(j, n_sides) + 1;
    faces = [faces; center_top, n_sides + j, n_sides + j_next];
end

end

function [vertices, faces] = create_box(p_min, p_max)
% Create a box mesh
x1 = p_min(1); y1 = p_min(2); z1 = p_min(3);
x2 = p_max(1); y2 = p_max(2); z2 = p_max(3);

vertices = [
    x1, y1, z1;  % 1
    x2, y1, z1;  % 2
    x2, y2, z1;  % 3
    x1, y2, z1;  % 4
    x1, y1, z2;  % 5
    x2, y1, z2;  % 6
    x2, y2, z2;  % 7
    x1, y2, z2;  % 8
];

faces = [
    % Bottom (y = y1)
    1, 2, 5;
    2, 6, 5;
    % Top (y = y2)
    3, 4, 7;
    4, 8, 7;
    % Front (z = z1)
    1, 4, 2;
    2, 4, 3;
    % Back (z = z2)
    5, 6, 8;
    6, 7, 8;
    % Left (x = x1)
    1, 5, 4;
    4, 5, 8;
    % Right (x = x2)
    2, 3, 6;
    3, 7, 6;
];

end

function write_stl_binary(filename, vertices, faces, solid_name)
% Write binary STL file
% vertices: Nx3 array of vertex coordinates
% faces: Mx3 array of face indices (1-based)
% solid_name: name for the solid (used in ASCII, ignored in binary but kept for reference)

if nargin < 4
    solid_name = 'matlab_export';
end

n_faces = size(faces, 1);

% Open file for binary writing
fid = fopen(filename, 'wb');
if fid == -1
    error('Could not open file for writing: %s', filename);
end

% Write 80-byte header
header = sprintf('%-80s', ['MATLAB STL: ' solid_name]);
fwrite(fid, header(1:80), 'char');

% Write number of triangles
fwrite(fid, n_faces, 'uint32');

% Write each triangle
for i = 1:n_faces
    % Get vertex indices
    v1 = vertices(faces(i,1), :);
    v2 = vertices(faces(i,2), :);
    v3 = vertices(faces(i,3), :);
    
    % Calculate normal
    edge1 = v2 - v1;
    edge2 = v3 - v1;
    normal = cross(edge1, edge2);
    norm_length = norm(normal);
    if norm_length > eps
        normal = normal / norm_length;
    else
        normal = [0, 0, 1];
    end
    
    % Write normal (3 floats)
    fwrite(fid, normal, 'float32');
    
    % Write vertices (3 x 3 floats)
    fwrite(fid, v1, 'float32');
    fwrite(fid, v2, 'float32');
    fwrite(fid, v3, 'float32');
    
    % Write attribute byte count (unused, set to 0)
    fwrite(fid, 0, 'uint16');
end

fclose(fid);

end