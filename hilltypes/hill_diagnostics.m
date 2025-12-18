close all; clc; clear;
hill_diagnostics_osimCurves

function hill_diagnostics_osimCurves()
% Hill-type diagnostics using each muscle's true OpenSim curves when available
% (Millard2012EquilibriumMuscle). Falls back to analytic curves if needed.

%% ===== USER INPUTS =====
modelPath = '/Users/chen/Downloads/MouseArmProject_v1_0/model/scaled_mouse.osim';
stoLen = '/Users/chen/Downloads/MouseArmProject_v1_0/analysis_and_processing_scripts/sim scripts/data/Predicted EMG/ScaleSet_MuscleAnalysis_FiberLength.sto';
stoVel = '/Users/chen/Downloads/MouseArmProject_v1_0/analysis_and_processing_scripts/sim scripts/data/Predicted EMG/ScaleSet_MuscleAnalysis_FiberVelocity.sto';
stoAct = '/Users/chen/Downloads/MouseArmProject_v1_0/model/subject01-scaled_StaticOptimization_activation.sto';

muscles = {'Biceps_Long_Head','Biceps_Short_Head', ...
           'Brachialis_Proximal_Head','Brachialis_Distal_Head', ...
           'Triceps_Long_Head','Triceps_Lat_Head','Anconeus'};

%% ===== Load model & per-muscle parameters/curves =====
import org.opensim.modeling.*
model = Model(modelPath); state = model.initSystem();
ms = model.getMuscles();

% storage
params = struct();   % one struct per muscle name

for i=0:ms.getSize()-1
    m  = ms.get(i);
    nm = char(m.getName());

    S.name   = nm;
    S.Fmax   = tryget(@() m.getMaxIsometricForce(), 1);     % N
    S.l0     = tryget(@() m.getOptimalFiberLength(), 0.02); % m
    S.vmax   = tryget(@() m.getMaxContractionVelocity(),10);% in l0/s
    S.phi0   = tryget(@() m.getPennationAngleAtOptimalFiberLength(), 0.0);

    % default (for non-Millard): analytic curves
    S.hasMillard = false;
    S.ActiveFL   = @(ln) max(0, 1 - ((ln-1)/0.25).^2);
    S.PassiveFL  = @(ln) (ln>1).* (exp(5*(ln-1))-1) * 0.03; % gentle exp
    S.FV         = @(vn) (vn>=0).*(1.8) + (vn<0).*((1-0.3*(-vn))./(1+1.5*(-vn)));

    % Try to cast to Millard and read true curves
    try
        mm = Millard2012EquilibriumMuscle.safeDownCast(m);
        if ~mm.isNull()
            S.hasMillard = true;
            % these curve objects expose calcValue(x)
            actFL  = mm.getActiveForceLengthCurve();
            pasFL  = mm.getFiberForceLengthCurve();      % passive fiber curve
            fv     = mm.getForceVelocityCurve();

            % Wrap into MATLAB function handles
            S.ActiveFL  = @(ln) max(0, actFL.calcValue(ln));   % normalized
            S.PassiveFL = @(ln) max(0, pasFL.calcValue(ln));   % normalized
            % OpenSim FV uses v_norm = v_f / (vmax*l0); shortening usually negative
            S.FV        = @(vn) max(0, fv.calcValue(vn));
        end
    catch
        % keep analytic fallback
    end

    params.(nm) = S;
end

%% ===== Read STOs; align activation time base =====
Tlen = readtable(stoLen,'FileType','text'); Tlen(:,all(ismissing(Tlen)))=[];
Tvel = readtable(stoVel,'FileType','text'); Tvel(:,all(ismissing(Tvel)))=[];
Tact = readtable(stoAct,'FileType','text'); Tact(:,all(ismissing(Tact)))=[];

t = Tlen{:,1};
if any(abs(Tvel{:,1}-t) > 1e-10)
    error('FiberLength and FiberVelocity times differ; resample first.');
end
if ~isequal(height(Tact),height(Tlen)) || any(abs(Tact{:,1}-t)>1e-10)
    % resample activation to length/velocity time base
    Tact_rs = table(t,'VariableNames',{'time'});
    for k=1:numel(muscles)
        m = muscles{k};
        Tact_rs.(m) = interp1(Tact{:,1}, Tact.(m), t, 'linear','extrap');
    end
    Tact = Tact_rs;
end

%% ===== Helper: pennation update φ(lf) (constant-thickness) =====
pennation = @(lf,l0,phi0) asin( min( max( (sin(phi0).*l0)./max(lf,1e-6), -1), 1) );

%% ===== Compute & plot =====
nm = numel(muscles);
figure('Name','Hill (OpenSim curves + activation + passive)','Color','w');

TotalPower = zeros(length(t),1);

for k=1:nm
    m  = muscles{k};
    P  = params.(m);

    % sanity on columns
    mustHave(Tlen,m); mustHave(Tvel,m); mustHave(Tact,m);

    lf = Tlen.(m);                  % m
    lv = Tvel.(m);                  % m/s
    a  = max(0,min(1,Tact.(m)));    % [0..1]

    ln = lf ./ P.l0;                % normalized length
    vn = lv ./ (P.vmax*P.l0);       % normalized velocity (OpenSim sign)

    % evaluate curves
    fl_act = P.ActiveFL(ln);        % normalized
    fl_pas = P.PassiveFL(ln);       % normalized (passive)
    fv     = P.FV(vn);              % normalized

    % fiber force (normalized + absolute)
    Fnorm_fiber = a .* fl_act .* fv + fl_pas;     % total normalized fiber force
    phi         = pennation(lf, P.l0, P.phi0);
    Fabs_mtu    = P.Fmax .* Fnorm_fiber .* cos(phi);  % project to tendon/MTU

    % mechanical power (approx; ignoring series compliance)
    Pm = Fabs_mtu .* lv;            % W
    TotalPower = TotalPower + Pm;

    % classify concentric/eccentric by fiber shortening/lengthening
    isShort = vn < 0;   % Millard uses vn<0 = shortening (concentric)
    isLen   = vn >= 0;

    % ---- F–L (active+passive shown) ----
    subplot(3,nm,k);
    scatter(ln, a.*fl_act.*fv, 10, a, 'filled'); hold on
    scatter(ln, fl_pas, 10, 'r', 'filled');
    grid on; xlabel('l_f / l_0'); ylabel('F/F_{max}');
    title(strrep([m '  (F–L)'],'_','\_'));
    legend('active (a·FL·FV)','passive FL','Location','best'); legend boxoff

    % ---- F–V (active only, colored by a; show ecc/conc markers) ----
    subplot(3,nm,k+nm);
    scatter(vn(isShort), a(isShort).*fl_act(isShort).*fv(isShort), 10, a(isShort),'filled'); hold on
    scatter(vn(isLen),   a(isLen).*fl_act(isLen).*fv(isLen),     10, a(isLen),'filled');
    grid on; xlabel('v_f / v_{max}'); ylabel('F/F_{max}');
    title(strrep([m '  (F–V, active only)'],'_','\_'));
    cb=colorbar; cb.Label.String='activation';
    yl = ylim; plot([0 0],yl,'k:'); ylim(yl);

    % ---- Absolute power over time ----
    subplot(3,nm,k+2*nm);
    plot(t, Pm, 'k-','LineWidth',1.0); grid on
    xlabel('time (s)'); ylabel('W');
    title(strrep([m '  power (fiber→MTU)'],'_','\_'));
end
sgtitle('Hill-type: per-muscle OpenSim curves (active+passive)');

figure('Color','w');
plot(t, TotalPower, 'LineWidth',1.4); grid on
xlabel('time (s)'); ylabel('W');
title('Total mechanical power (selected muscles)');
end

%% ===== utilities =====
function x = tryget(fun, fallback)
try, x = fun(); catch, x = fallback; end
end
function mustHave(T, nm)
assert(ismember(nm, T.Properties.VariableNames), "Column not found: %s", nm);
end