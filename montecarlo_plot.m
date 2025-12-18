files = dir('synthetic6dof_*.mot');

allW_conc = [];
allW_ecc  = [];

for i=1:length(files)
    tbl = readtable(files(i).name,'FileType','text','Delimiter','\t','HeaderLines',7);
    tbl.Properties.VariableNames = {'time','elv_angle','extension_angle','rotation_angle','elbow_flex','wrist_angle','radius_rot'};

    [W_conc, W_ecc] = estimateTorqueWork(tbl); % where this is your torque function
    allW_conc(end+1) = W_conc;
    allW_ecc(end+1)  = W_ecc;
end

figure; 
histogram(allW_conc,20); hold on; histogram(allW_ecc,20);
legend('Concentric Work','Eccentric Work');
title('Monte Carlo Work Distribution across motor repertoires');
xlabel('Joules');