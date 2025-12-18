% coordinate names (same order as model)
coords = {'elv_angle','extension_angle','rotation_angle', ...
          'elbow_flex','wrist_angle','radius_rot'};

% time vector
T = linspace(0,2,501)'; % 2 seconds, 500 samples

% multi sweep joint trajectories
elv_angle       = 20*sin(2*pi*0.5*T);  % deg
extension_angle = 10*sin(2*pi*0.7*T+0.5);
rotation_angle  = 15*sin(2*pi*0.4*T+1.0);
elbow_flex      = 45*sin(2*pi*0.6*T+1.5);
wrist_angle     = 30*sin(2*pi*0.9*T+0.2);
radius_rot      = 25*cos(2*pi*0.3*T); 

data = [T elv_angle extension_angle rotation_angle elbow_flex wrist_angle radius_rot];

% --- write MOT file --- %
fname='synthetic6dof.mot';
fid=fopen(fname,'w');

fprintf(fid,'Synthetic Multi Sweep Motion\n');
fprintf(fid,'version=1\n');
fprintf(fid,'nRows=%d\n',size(data,1));
fprintf(fid,'nColumns=7\n');
fprintf(fid,'inDegrees=yes\n');
fprintf(fid,'endheader\n');

% write header line for columns
fprintf(fid,'time\t%s\t%s\t%s\t%s\t%s\t%s\n',coords{:});

fclose(fid);
writematrix(data,fname,'Delimiter','\t','WriteMode','append','FileType','text');
disp(['written: ' fname]);