% number of samples
N = 200; 

% coordinate names
coords = {'elv_angle','extension_angle','rotation_angle', ...
          'elbow_flex','wrist_angle','radius_rot'};

for k = 1:N

    % random parameters each run
    freq = 0.4 + 0.8*rand(1,6);        % random frequencies 0.4-1.2 Hz
    amp  = [20 10 15 45 30 25] .* (0.8 + 0.4*rand(1,6));   % amplitude jitter ±20%
    phi  = 2*pi*rand(1,6);             % random phase offsets

    % time
    T = linspace(0,2,501)';

    % synth motion
    X = zeros(length(T),6);
    for j=1:6
        X(:,j) = amp(j)*sin(2*pi*freq(j)*T + phi(j));
    end

    data = [T X];

    fname = sprintf('synthetic6dof_%03d.mot',k);

    fid=fopen(fname,'w');
    fprintf(fid,'Synthetic Multi Sweep Motion %d\n',k);
    fprintf(fid,'version=1\n');
    fprintf(fid,'nRows=%d\n',size(data,1));
    fprintf(fid,'nColumns=7\n');
    fprintf(fid,'inDegrees=yes\n');
    fprintf(fid,'endheader\n');
    fprintf(fid,'time\t%s\t%s\t%s\t%s\t%s\t%s\n',coords{:});
    fclose(fid);

    writematrix(data,fname,'Delimiter','\t','WriteMode','append','FileType','text');

end

disp('DONE generating MOTs');