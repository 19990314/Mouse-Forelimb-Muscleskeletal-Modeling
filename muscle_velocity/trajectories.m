import org.opensim.modeling.*

%% === Workspace cloud from many MOT files (marker: 'hand') ===
% CONFIG
modelPath   = '/Users/chen/Downloads/MouseArmProject_v1_0/model/scaled_mouse.osim';   % <-- your scaled model
motPattern  = 'synthetic6dof_*.mot';       % files to process
markerName  = 'handm';                      % exact marker name
coords      = {'elv_angle','extension_angle','rotation_angle', ...
               'elbow_flex','wrist_angle','radius_rot'};    % model coord names
inDegrees   = true;                        % your .mot header: inDegrees=yes

% ---- OpenSim import
try
    import org.opensim.modeling.*
catch
    error(['OpenSim MATLAB API not found. In OpenSim GUI: Help > Scripting with MATLAB, ', ...
           'or add opensim-modeling.jar to your Java classpath.']);
end

% ---- Load model
model = Model(modelPath);
state = model.initSystem();


% Handle to marker
try
    mkr = model.getMarkerSet().get(markerName);
catch
    error('Marker "%s" not found in model. Check exact spelling.', markerName);
end

% Quick coord existence check
cs = model.getCoordinateSet();
for k = 1:numel(coords)
    if ~cs.contains(coords{k})
        error('Model has no coordinate named "%s".', coords{k});
    end
end

% ---- Collect point cloud
files = dir(motPattern);
XYZall = [];                 % big [N x 3] cloud
fileId = [];                 % optional: which file each point came from

for f = 1:numel(files)
    fn = files(f).name;
    % Read MOT (skip 7-line header typical for OpenSim)
    T = readtable(fn, 'FileType','text','Delimiter','\t','HeaderLines',7);

    % If readtable didn't capture names, set them now (time + 6 coords)
    if width(T) >= 7
        T.Properties.VariableNames(1:7) = {'time', coords{:}};
    else
        error('File %s has fewer than 7 columns.', fn);
    end

    % Convert to radians if needed
    Q = zeros(height(T), numel(coords));
    for k = 1:numel(coords)
        ang = T.(coords{k});
        if inDegrees
            Q(:,k) = deg2rad(ang);
        else
            Q(:,k) = ang;
        end
    end
        

    % For each row: set coordinates -> realize -> marker position
    for i = 1:height(T)
        % set each coordinate value
        for k = 0:cs.getSize()-1
            % if not in our list, leave as current; set only those we have
        end
        for k = 1:numel(coords)
            cs.get(coords{k}).setValue(state, Q(i,k), false);
        end
        % realize to Position stage
        model.realizePosition(state);

        % marker position in ground
        p = mkr.getLocationInGround(state);  % Vec3
        XYZall(end+1, :) = [p.get(0), p.get(1), p.get(2)]; %#ok<AGROW>
        fileId(end+1, 1) = f;                 %#ok<AGROW>
    end

    fprintf('Processed %s (%d frames)\n', fn, height(T));
end

% ---- Plot workspace cloud
figure; 
scatter3(XYZall(:,1), XYZall(:,2), XYZall(:,3), 6, fileId, 'filled'); 
axis equal; grid on
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Hand marker workspace cloud (%d motions, %d points)', numel(files), size(XYZall,1)));
colorbar; view(135,20);

% Optional: density via alpha shape (quick hull volume estimate)
try
    shp = alphaShape(XYZall, inf);
    V = volume(shp);
    fprintf('Approx workspace volume = %.3e m^3\n', V);
catch
    % alphaShape requires MATLAB's geometry toolbox; safe to skip
end