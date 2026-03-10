%% PS3 Problem 2a: f-I curve (no refractory), theory + simulation
clear; clc;

Rm    = 100e6;
Cm    = 100e-12;
Erest = -70e-3;
Vth   = -60e-3;
Vreset = Erest;
tau = Rm*Cm;
dt = 0.1e-3;
T  = 0.5;  
t  = 0:dt:T;

% Theory: 0..1000 pA in 1 pA steps
I_theory = (0:1000)*1e-12;           % A
Vinf = Erest + Rm*I_theory;          % V
f_theory = zeros(size(I_theory));

idx = Vinf > Vth;
f_theory(idx) = 1 ./ ( tau .* log( (Vreset - Vinf(idx)) ./ (Vth - Vinf(idx)) ) );

figure; hold on;
plot(I_theory*1e12, f_theory, 'LineWidth', 1.6);

% Simulation: 0..1000 pA in 50 pA steps
I_sim = (0:50:1000)*1e-12;
f_sim = zeros(size(I_sim));

for k=1:numel(I_sim)
    Iext = I_sim(k)*ones(size(t));
    [~, spk] = localSimSlide(t, dt, Rm, Cm, Erest, Vth, Vreset, Iext);
    f_sim(k) = spk / T;
end

plot(I_sim*1e12, f_sim, 'o', 'LineWidth', 1.6);

xlabel('I_{ext} (pA)'); ylabel('Firing rate (Hz)');
title('f–I curve (no refractory): theory + simulation');
legend('Theory (1 pA steps)','Simulation (50 pA steps)','Location','best');
grid on;

function [Vm, nSpikes] = localSimSlide(t, dt, Rm, Cm, Erest, Vth, Vreset, Iext)
    Vm = zeros(size(t));
    Vm(1) = Vreset;
    nSpikes = 0;

    for it = 1:length(t)-1
        Vm(it+1) = Vm(it) + dt*(Iext(it)*Rm + Erest - Vm(it)) / (Rm*Cm);
        if Vm(it+1) >= Vth
            nSpikes = nSpikes + 1;
            Vm(it+1) = Vreset;
        end
    end
end