%% PS3 Problem 2b: f-I curve WITH absolute refractory (Tref = 2 ms)
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

% Absolute refractory
Tref = 2e-3;                  % 2 ms
Nref = ceil(Tref/dt);         % number of time steps to hold in refractory

%% ---- Theory: 0..1000 pA step 1 pA ----
I_theory = (0:1000)*1e-12;             % A
Vinf = Erest + Rm*I_theory;            % V

f_theory = zeros(size(I_theory));
idx = Vinf > Vth;                       % spiking condition

Tisi = tau .* log( (Vreset - Vinf(idx)) ./ (Vth - Vinf(idx)) );  % from Problem 1(a)
f_theory(idx) = 1 ./ (Tisi + Tref);     % refractory adds extra dead time

figure; hold on;
plot(I_theory*1e12, f_theory, 'LineWidth', 1.8);

%% ---- Simulation: 0..1000 pA step 50 pA ----
I_sim = (0:50:1000)*1e-12;             % A
f_sim = zeros(size(I_sim));

for k = 1:numel(I_sim)
    Iext = I_sim(k)*ones(size(t));     % constant current
    [~, spk] = localSimSlideRef(t, dt, Rm, Cm, Erest, Vth, Vreset, Iext, Nref);
    f_sim(k) = spk / T;
end

plot(I_sim*1e12, f_sim, 'o', 'LineWidth', 1.8);

xlabel('I_{ext} (pA)');
ylabel('Firing rate (Hz)');
title('f–I curve (T_{ref}=2 ms): theory + simulation');
legend('Theory (1 pA steps)', 'Simulation (50 pA steps)', 'Location', 'best');
grid on;

%% ---- Local function: slide update + threshold/reset + refractory counter ----
function [Vm, nSpikes] = localSimSlideRef(t, dt, Rm, Cm, Erest, Vth, Vreset, Iext, Nref)
    Vm = zeros(size(t));
    Vm(1) = Vreset;
    nSpikes = 0;

    ref_count = 0;

    for it = 1:length(t)-1
        if ref_count > 0
            Vm(it+1) = Vreset;
            ref_count = ref_count - 1;
        else
            Vm(it+1) = Vm(it) + dt*(Iext(it)*Rm + Erest - Vm(it)) / (Rm*Cm);

            if Vm(it+1) >= Vth
                nSpikes = nSpikes + 1;
                Vm(it+1) = Vreset;
                ref_count = Nref;
            end
        end
    end
end