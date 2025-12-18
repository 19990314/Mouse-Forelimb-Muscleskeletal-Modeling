function [W_conc, W_ecc] = estimateTorqueWork(tbl)

t   = tbl.time(:);
qdeg= tbl.elbow_flex(:);
q    = deg2rad(qdeg);
dt   = [diff(t); median(diff(t))];
qd   = [diff(q)./diff(t); 0];

% model params
r_flex=0.0020; r_ext=0.0018;
Fmax_flex=0.8; Fmax_ext=0.9;
L0_flex=0.008; L0_ext=0.008; 
vmax_flex=10; vmax_ext=10;
q0_flex=deg2rad(40); q0_ext=deg2rad(40);

% lengths & vel
dL_flex=-r_flex*(q-q0_flex); dL_ext=r_ext*(q-q0_ext);
Lmt_flex=L0_flex+dL_flex; Lmt_ext=L0_ext+dL_ext;
l_flex=Lmt_flex/L0_flex; l_ext=Lmt_ext/L0_ext;
ldot_f=-r_flex*qd; ldot_e=r_ext*qd;
v_flex=ldot_f/(vmax_flex*L0_flex); v_ext=ldot_e/(vmax_ext*L0_ext);

% curves
w=0.25; fl=@(l)exp(-((l-1)./w).^2);
fv=@(v)max(0.1, min(1.8,(1-0.3*v)./(1+1.5*v)));
fpass=@(l)(l>1).*(exp(5*(l-1))-1)*0.02;

% activation
a_flex=0.3+0.4*(qd>0)+0.1*(qd==0);
a_ext =0.3+0.4*(qd<0)+0.1*(qd==0);

% forces & torque
F_flex=a_flex.*Fmax_flex.*fl(l_flex).*fv(v_flex)+Fmax_flex.*fpass(l_flex);
F_ext =a_ext .*Fmax_ext .*fl(l_ext ).*fv(v_ext )+Fmax_ext .*fpass(l_ext);
tau=r_flex.*F_flex - r_ext.*F_ext;

ecc  = qd<0;
conc = qd>0;

dtheta=[diff(q);0];
inst_work=tau.*dtheta;

W_conc=sum(inst_work(conc),'omitnan');
W_ecc =sum(inst_work(ecc ),'omitnan');
end