% load mot
tbl = readtable('synthetic6dof.mot','FileType','text','Delimiter','\t','HeaderLines',7);
tbl.Properties.VariableNames = { ...
    'time', ...
    'elv_angle', ...
    'extension_angle', ...
    'rotation_angle', ...
    'elbow_flex', ...
    'wrist_angle', ...
    'radius_rot' };

% sample moment arm atlas (toy demonstration: radial workspace)
% workspace projection (elbow_flex vs wrist_angle)
figure; 
scatter(tbl.elbow_flex, tbl.wrist_angle,10,T,'filled');
xlabel('Elbow flex (deg)'); ylabel('Wrist angle (deg)');
title('Forelimb Workspace Map');
colorbar; colormap turbo;

% crude eccentric vs concentric segmentation
eccentric = tbl.elbow_flex < 0;
concentric = tbl.elbow_flex >=0;

figure;
plot(tbl.time,tbl.elbow_flex,'k','LineWidth',1.5); hold on;
plot(tbl.time(eccentric), tbl.elbow_flex(eccentric),'r.','MarkerSize',12)
plot(tbl.time(concentric),tbl.elbow_flex(concentric),'b.','MarkerSize',12)

legend('angle','eccentric phase','concentric phase','Location','best')
xlabel('time (s)'); ylabel('degrees');
title('Elbow Flexion Cycle: Concentric vs Eccentric Classification')
grid on