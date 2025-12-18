% --------- USER PATHS ---------
stoForce = '/Users/chen/Downloads/MouseArmProject_v1_0/analysis_and_processing_scripts/sim scripts/data/Predicted EMG/ScaleSet_MuscleAnalysis_FiberForce.sto';
stoLen   = '/Users/chen/Downloads/MouseArmProject_v1_0/analysis_and_processing_scripts/sim scripts/data/Predicted EMG/ScaleSet_MuscleAnalysis_FiberLength.sto';
stoVel   = '/Users/chen/Downloads/MouseArmProject_v1_0/analysis_and_processing_scripts/sim scripts/data/Predicted EMG/ScaleSet_MuscleAnalysis_FiberVelocity.sto';

% muscles to analyze
muscles = {'Biceps_Short_Head','Triceps_Long_Head'};

% --------- LOAD STOs ----------
tblForce = readtable(stoForce,'FileType','text','Delimiter','\t','HeaderLines',11);
tblLen   = readtable(stoLen,  'FileType','text','Delimiter','\t','HeaderLines',11);
tblVel   = readtable(stoVel,  'FileType','text','Delimiter','\t','HeaderLines',11);

% --------- NORMALIZE ----------
for m = 1:length(muscles)
    mus = muscles{m};
    
    F = tblForce.(mus);
    L = tblLen.(mus);
    V = tblVel.(mus);
    
    Lnorm = L ./ mean(L);          % normalized fiber length
    Vnorm = V ./ max(abs(V));      % normalized fiber velocity
    
    figure;
    scatter(Lnorm,F,12,'filled')
    xlabel('Normalized Length')
    ylabel('Force (N)')
    title(['Force-Length: ' mus])
    grid on
    saveas(gcf,[mus '_Hill_FL.png'])
    
    figure;
    scatter(Vnorm,F,12,'filled')
    xlabel('Normalized Velocity')
    ylabel('Force (N)')
    title(['Force-Velocity: ' mus])
    grid on
    saveas(gcf,[mus '_Hill_FV.png'])
end