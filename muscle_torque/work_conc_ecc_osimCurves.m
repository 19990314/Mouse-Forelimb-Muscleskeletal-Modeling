function Tsum = work_conc_ecc_osimCurves()
% Concentric vs Eccentric work per muscle using OpenSim (Millard) curves.
% Outputs a table and figures; also writes 'muscle_work_summary.csv'.

%% ===== USER INPUTS =====
modelPath = '/Users/chen/Downloads/MouseArmProject_v1_0/model/scaled_mouse.osim';
stoLen = '/Users/chen/Downloads/MouseArmProject_v1_0/analysis_and_processing_scripts/sim scripts/data/Predicted EMG/ScaleSet_MuscleAnalysis_FiberLength.sto';
stoVel = '/Users/chen/Downloads/MouseArmProject_v1_0/analysis_and_processing_scripts/sim scripts/data/Predicted EMG/ScaleSet_MuscleAnalysis_FiberVelocity.sto';
stoAct = '/Users/chen/Downloads/MouseArmProject_v1_0/model/subject01-scaled_StaticOptimization_activation.sto';

muscles = {'Pectoralis_Major_Posterior','Pectoralis_Major_Anterior', ...
           'Latissimus_Dorsi_Caudal','Latissimus_Dorsi_Rostral', ...
           'Triceps_Long_Head','Biceps_Short_Head','Deltoid_Medial_CHK_MMTARM'};

%% ===== Load model & muscle parameters/curves =====
import org.opensim.modeling.*
model = Model(modelPath); state = model.initSystem();
ms = model.getMuscles();

P = struct();   % per-muscle param & curve handles
for i=0:ms.getSize()-1
    m  = ms.get(i); nm = char(m.getName());
    S.name = nm;
    S.Fmax = tryget(@() m.getMaxIsometricForce(), 1);       % N
    S.l0   = tryget(@() m.getOptimalFiberLength(), 0.02);   % m
    S.vmax = tryget(@() m.getMaxContractionVelocity(), 10); % in l0/s
    S.phi0 = tryget(@() m.getPennationAngleAtOptimalFiberLength(), 0.0);

    % defaults (if not Millard)
    S.ActiveFL  = @(ln) max(0, 1 - ((ln-1)/0.25).^2);
    S.PassiveFL = @(ln) (ln>1).* (exp(5*(ln-1))-1) * 0.03;
    S.FV        = @(vn) (vn>=0).*(1.8) + (vn<0).*((1-0.3*(-vn))./(1+1.5*(-vn)));

    % try Millard
    try
        mm = Millard2012EquilibriumMuscle.safeDownCast(m);
        if ~mm.isNull()
            actFL = mm.getActiveForceLengthCurve();
            pasFL = mm.getFiberForceLengthCurve();
            fv    = mm.getForceVelocityCurve();
            S.ActiveFL  = @(ln) max(0, actFL.calcValue(ln));
            S.PassiveFL = @(ln) max(0, pasFL.calcValue(ln));
            S.FV        = @(vn) max(0, fv.calcValue(vn));
        end
    catch
        % keep defaults
    end

    P.(nm) = S;
end

%% ===== Read & align STO tables =====
Tlen = readtable(stoLen,'FileType','text'); Tlen(:,all(ismissing(Tlen)))=[];
Tvel = readtable(stoVel,'FileType','text'); Tvel(:,all(ismissing(Tvel)))=[];
Tact = readtable(stoAct,'FileType','text'); Tact(:,all(ismissing(Tact)))=[];
t = Tlen{:,1};
if any(abs(Tvel{:,1}-t) > 1e-10)
    error('FiberLength and FiberVelocity have different time bases.');
end
if ~isequal(height(Tact),height(Tlen)) || any(abs(Tact{:,1}-t)>1e-10)
    Tact_rs = table(t,'VariableNames',{'time'});
    for k=1:numel(muscles)
        m = muscles{k};
        assert(ismember(m,Tact.Properties.VariableNames),"Missing activation for %s",m);
        Tact_rs.(m) = interp1(Tact{:,1}, Tact.(m), t, 'linear','extrap');
    end
    Tact = Tact_rs;
end

%% ===== Helpers =====
pennation = @(lf,l0,phi0) asin( min( max( (sin(phi0).*l0)./max(lf,1e-6), -1), 1) );

%% ===== Loop muscles: power & work =====
nm = numel(muscles);
Wconc      = zeros(nm,1);  % including passive
Wecc       = zeros(nm,1);
Wconc_act  = zeros(nm,1);  % active-only
Wecc_act   = zeros(nm,1);

for k=1:nm
    m  = muscles{k};
    assert(ismember(m,Tlen.Properties.VariableNames),"No length for %s",m);
    assert(ismember(m,Tvel.Properties.VariableNames),"No velocity for %s",m);
    assert(ismember(m,Tact.Properties.VariableNames),"No activation for %s",m);

    pk = P.(m);

    lf = Tlen.(m);              % m
    lv = Tvel.(m);              % m/s  (Millard: vn<0 = shortening)
    a  = max(0,min(1,Tact.(m)));

    ln = lf ./ pk.l0;
    vn = lv ./ (pk.vmax*pk.l0);

    fl_act = pk.ActiveFL(ln);
    fl_pas = pk.PassiveFL(ln);
    fv     = pk.FV(vn);

    Fnorm_active = a .* fl_act .* fv;      % normalized
    Fnorm_total  = Fnorm_active + fl_pas;  % + passive
    phi          = pennation(lf, pk.l0, pk.phi0);

    Fmtu_active = pk.Fmax .* Fnorm_active .* cos(phi); % N
    Fmtu_total  = pk.Fmax .* Fnorm_total  .* cos(phi); % N

    % power convention:
    % concentric power (delivered by muscle) >=0 when vn<0 (shortening)
    Pconc_total = -Fmtu_total .* lv .* (vn<0);
    Pecc_total  =  Fmtu_total .* lv .* (vn>=0);
    Pconc_act   = -Fmtu_active.* lv .* (vn<0);
    Pecc_act    =  Fmtu_active.* lv .* (vn>=0);

    % integrate with trapz (Joules)
    Wconc(k)     = trapz(t, Pconc_total);
    Wecc(k)      = trapz(t, Pecc_total);
    Wconc_act(k) = trapz(t, Pconc_act);
    Wecc_act(k)  = trapz(t, Pecc_act);
end

%% ===== Summaries & plots =====
Tsum = table(muscles(:), Wconc, Wecc, Wconc_act, Wecc_act, ...
    'VariableNames', {'muscle','W_conc_total','W_ecc_total','W_conc_active','W_ecc_active'});

writetable(Tsum,'muscle_work_summary.csv');

% Bar plots (total)
figure('Color','w');
subplot(1,2,1);
barh(Wconc); set(gca,'YTick',1:nm,'YTickLabel',strrep(muscles,'_','\_'));
xlabel('J'); title('Concentric work (total)'); grid on
subplot(1,2,2);
barh(Wecc); set(gca,'YTick',1:nm,'YTickLabel',strrep(muscles,'_','\_'));
xlabel('J'); title('Eccentric work (total)'); grid on

% Stacked total across muscles
figure('Color','w');
totC = sum(Wconc); totE = sum(Wecc);
bar([totC totE]); set(gca,'XTickLabel',{'Concentric','Eccentric'});
ylabel('J'); title('Total work (all selected muscles)'); grid on

% Active-only stacked
figure('Color','w');
totC_a = sum(Wconc_act); totE_a = sum(Wecc_act);
bar([totC_a totE_a]); set(gca,'XTickLabel',{'Concentric (active)','Eccentric (active)'});
ylabel('J'); title('Active-only work'); grid on

disp('Wrote: muscle_work_summary.csv');
end

%% ===== utilities =====
function x = tryget(fun, fallback)
try
    x = fun(); catch, x = fallback; end
end



%% ===== fib =====


muscleGroups = {
    'Shoulder', ...
    'Elbow', ...
    'Elbow', ...
    'Shoulder', ...
    'Elbow', ...
    'Elbow', ...
    'Elbow', ...
    'Elbow', ...
    'Elbow', ...
    'Shoulder', ...
    'Elbow', ...
    'Shoulder', ...
    'Shoulder', ...
    'Wrist', ...
    'Wrist', ...
    'Wrist', ...
    'Elbow', ...
    'Shoulder', ...
    'Shoulder', ...
    'Shoulder', ...
    'Shoulder'
};
figure('Position', [100 100 700 800]);

% Heatmap
subplot(4,1,1:2);
imagesc(fiberVelocities, [-5, 5]);
colormap(parula);
colorbar('SouthOutside', 'Ticks', -5:1:5, 'TickLabels', -5:1:5);
xlabel('Postures'); title('Fiber Velocity (fiber lengths/s)');

% Label Y-axis with color-coded groups
yticks(1:length(muscleNames));
yticklabels(muscleNames);
set(gca, 'YTickLabelRotation', 0);

% Color group labels
for i = 1:length(muscleNames)
    if strcmp(group{i}, 'Shoulder')
        text(-2, i, muscleNames{i}, 'Color', [0.2 0.2 1], 'HorizontalAlignment', 'right');
    elseif strcmp(group{i}, 'Elbow')
        text(-2, i, muscleNames{i}, 'Color', [1 0.2 0.2], 'HorizontalAlignment', 'right');
    elseif strcmp(group{i}, 'Wrist')
        text(-2, i, muscleNames{i}, 'Color', [0.1 0.8 0.1], 'HorizontalAlignment', 'right');
    end
end

% Add horizontal lines to separate groups
hold on;
elbowIdx = find(strcmp(group, 'Elbow'), 1);
wristIdx = find(strcmp(group, 'Wrist'), 1);
yline(elbowIdx - 0.5, '--k');
yline(wristIdx - 0.5, '--k');


% Function to extract indices and plot
muscleGroupNames = {'Shoulder', 'Elbow', 'Wrist'};
colors = {[0.2 0.2 1], [1 0.2 0.2], [0.1 0.8 0.1]};

for g = 1:3
    subplot(4,1,2+g);
    idx = strcmp(group, muscleGroupNames{g});
    data = fiberVelocities(idx, :);
    
    meanV = mean(data,1);
    stdV = std(data,[],1);
    
    fill([1:nPostures, fliplr(1:nPostures)], ...
         [meanV+stdV, fliplr(meanV-stdV)], ...
         colors{g}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on;
    plot(meanV, 'Color', colors{g}, 'LineWidth', 2);
    
    ylim([-15, 15]);
    yline(0, '--');
    title(sprintf('%s Muscles', muscleGroupNames{g}));
    ylabel('Fiber Velocity (l/s)');
end
xlabel('Postures');