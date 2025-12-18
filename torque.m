%% === Inputs: your kinematics table (angles in degrees) ===
% If you need column names, set them explicitly (uncomment next block)
% tbl = readtable('synthetic6dof.mot','FileType','text','Delimiter','\t','HeaderLines',7);
% tbl.Properties.VariableNames = {'time','elv_angle','extension_angle','rotation_angle', ...
%                                 'elbow_flex','wrist_angle','radius_rot'};

t   = tbl.time(:);
qdeg= tbl.elbow_flex(:);              % elbow flexion angle (deg)
q    = deg2rad(qdeg);                 % radians
dt   = [diff(t); median(diff(t))];    % robust dt
qd   = [diff(q)./diff(t); 0];         % angular velocity (rad/s)

%% === Tiny Hill-type pair for the elbow (flexor & extensor) ===
% Moment arms (meters) — rough mouse-scale; adjust with your model's values
r_flex = 0.0020;   % biceps-like elbow flexor
r_ext  = 0.0018;   % triceps-like elbow extensor

% Muscle architectural params (per muscle)
Fmax_flex = 0.8;   % max isometric force [N] (toy mouse-scale)
Fmax_ext  = 0.9;   % [N]
L0_flex   = 0.008; % optimal fiber length [m]
L0_ext    = 0.008; % [m]
vmax_flex = 10;    % vmax ~ 10 * L0 / s (normalized)
vmax_ext  = 10;

% Reference angle where fibers are near optimal length (rad)
q0_flex = deg2rad(40);  % flexor optimal joint angle
q0_ext  = deg2rad(40);  % extensor

% ----- Kinematics → fiber length & velocity (per muscle) -----
% For a single-DOF hinge, fiber length change ≈ -r * (q - q0) for the agonist
dL_flex   = -r_flex*(q - q0_flex);          % [m]
dL_ext    =  r_ext *(q - q0_ext );          % [m]
Lmt_flex  = L0_flex + dL_flex;              % crude MTU fiber length proxy
Lmt_ext   = L0_ext  + dL_ext;

ltil_flex = Lmt_flex./L0_flex;              % normalized length
ltil_ext  = Lmt_ext ./L0_ext;

% Fiber velocity: ldot ≈ -r*qdot (flexor), +r*qdot (extensor)
ldot_flex = -r_flex*qd;                     % [m/s]
ldot_ext  =  r_ext *qd;

vtil_flex = ldot_flex./(vmax_flex*L0_flex); % normalized fiber velocity
vtil_ext  = ldot_ext ./(vmax_ext *L0_ext);

% ----- Force–length (Gaussian) -----
w = 0.25; % width
fl = @(ltil) exp(-((ltil-1)./w).^2);

% ----- Force–velocity (smooth, clipped) -----
% Convention: vtil < 0 → shortening (concentric); vtil > 0 → lengthening (eccentric)
fv = @(vtil) max(0.1, min(1.8, (1 - 0.3*vtil)./(1 + 1.5*vtil)));  % simple Hill-like

% ----- Passive force (exponential beyond ~1.0 L0) -----
fpass = @(ltil) (ltil>1).* (exp(5*(ltil-1)) - 1)*0.02;

% ----- Phase-dependent activations (toy co-contraction rule) -----
% When qd>0 (flexing): more flexor activation; when qd<0: more extensor (braking)
a_flex = 0.3 + 0.4*(qd>0) + 0.1*(qd==0);    % 0.3 base + 0.4 during concentric flexion
a_ext  = 0.3 + 0.4*(qd<0) + 0.1*(qd==0);    % 0.3 base + 0.4 during eccentric (extension direction)

% ----- Muscle forces -----
F_flex = a_flex.*Fmax_flex .* fl(ltil_flex).*fv(vtil_flex) + Fmax_flex.*fpass(ltil_flex);
F_ext  = a_ext .*Fmax_ext  .* fl(ltil_ext ).*fv(vtil_ext ) + Fmax_ext .*fpass(ltil_ext );

% ----- Net joint torque (about elbow_flex) -----
tau = r_flex.*F_flex - r_ext.*F_ext;        % [N·m]

%% === Concentric vs eccentric labeling (w.r.t. elbow flexors) ===
% For flexors: concentric when qdot>0 (shortening), eccentric when qdot<0
eccentric   = qd < 0;
concentric  = qd > 0;

%% === Work (separate conc/ecc) ===
% Work = ∫ tau * dq
dtheta      = [diff(q); 0];
inst_work   = tau .* dtheta;                % incremental work [J] (N·m·rad)
% --- helpers (safe even if some NaNs exist) ---
omit = @(x) sum(x,'omitnan');     % fallback for nansum
sgn  = @(x) gradient(x);          % simple derivative

t   = tbl.time(:);
qdeg= tbl.elbow_flex(:);          % degrees
q    = deg2rad(qdeg);             
qd   = [diff(q)./diff(t); 0];     % rad/s

% parameters (toy defaults)
r_flex=0.0020; r_ext=0.0018; 
Fmax_flex=0.8; Fmax_ext=0.9; 
L0_flex=0.008; L0_ext=0.008; vmax_flex=10; vmax_ext=10;
q0_flex=deg2rad(40); q0_ext=deg2rad(40);

% lengths/velocities
dL_flex = -r_flex*(q - q0_flex);   dL_ext  =  r_ext*(q - q0_ext);
Lmt_flex= L0_flex + dL_flex;       Lmt_ext =  L0_ext + dL_ext;
l_flex  = Lmt_flex/L0_flex;        l_ext   =  Lmt_ext /L0_ext;
ldot_f  = -r_flex*qd;              ldot_e  =  r_ext*qd;
v_flex  = ldot_f /(vmax_flex*L0_flex);
v_ext   = ldot_e /(vmax_ext *L0_ext);

% curves
w=0.25; fl = @(l) exp(-((l-1)./w).^2);
fv = @(v) max(0.1, min(1.8, (1 - 0.3*v)./(1 + 1.5*v)));
fpass = @(l) (l>1).*(exp(5*(l-1)) - 1)*0.02;

% simple phase activations
a_flex = 0.3 + 0.4*(qd>0) + 0.1*(qd==0);
a_ext  = 0.3 + 0.4*(qd<0) + 0.1*(qd==0);

% forces & torque
F_flex = a_flex.*Fmax_flex.*fl(l_flex).*fv(v_flex) + Fmax_flex.*fpass(l_flex);
F_ext  = a_ext .*Fmax_ext .*fl(l_ext ).*fv(v_ext ) + Fmax_ext .*fpass(l_ext );
tau    = r_flex.*F_flex - r_ext.*F_ext;   % N·m

% concentric/eccentric labels for flexors
eccentric  = qd < 0;
concentric = qd > 0;

% work (use sum(...,'omitnan') instead of nansum)
dtheta    = [diff(q); 0];                 % rad
inst_work = tau .* dtheta;                % J (N·m·rad)

W_conc = omit(inst_work(concentric));
W_ecc  = omit(inst_work(eccentric));

fprintf('Elbow: Concentric work = %.4e J, Eccentric work = %.4e J\n', W_conc, W_ecc);

% quick plot
figure;
subplot(3,1,1); plot(t,qdeg,'k','LineWidth',1.4); grid on; ylabel('Elbow (deg)');
subplot(3,1,2); plot(t,tau,'LineWidth',1.6); grid on; ylabel('\tau (N·m)');
subplot(3,1,3); plot(t,qdeg,'k'); hold on
plot(t(eccentric),qdeg(eccentric),'r.',t(concentric),qdeg(concentric),'b.')
grid on; xlabel('Time (s)'); ylabel('deg'); legend('ang','ecc','conc'); 