function [timevector, v_1, v_2] = ML_coupl_2009(T, dt, Gsyn, Gh)

%%% Function that received maximal conductances for inhibitory synapse and
%%% h current and returns two voltage traces
%%% cell 1 is formerly biological, cell 2 is formerly model

% T should be in ms
% dt should be in ms
% Gsyn should be in uS
% Gh should be in uS

% parameters
C_m = 10.0;  % cell capacitance, nF
E_Ca = +100.0; % mV
tau_m_Ca = 0.001;  % ms
tau_0_K = 1.0/0.002; % ms
syn_tau = 100; % ms
s_r_1 = 7;  % mV
v_h_half = -50;  % mV
c_r= .00033; % 1/ms
v_kr = -110; % mV
e_h = -10; % mV
s_kr = -13; % mV
Vhalf_Ca = -20; % mV
Vslope_Ca = 15; % mV
Vhalf_K = -20; % mV
Vslope_K = 15; % mV
Gleak =200e-3; % uS
E_leak = -60;  % mV
Vhalf_syn = -25;  % mV
Vslope_syn = 15;  % mV
E_K = -80;  % mV
E_syn_12 = -90; % mV, from mod to bio, should  be -90
E_syn_21 = -70; % mV, from bio to mod should be -70
GCa = 200e-3;  % uS
GK = 200e-3;  % uS

% timeline definition
timevector = (0:dt:T-dt);
n_time = length(timevector);

% allocate stuff we want to save
v_1 = zeros(1,n_time);  % mV
v_2 = zeros(1,n_time);
% Ca
m_Ca_1 = zeros(1,n_time);  % pure
m_Ca_2 = zeros(1,n_time);
g_Ca_1 = zeros(1,n_time);
g_Ca_2 = zeros(1,n_time);
% K
m_K_1= zeros(1,n_time);  % pure
m_K_2= zeros(1,n_time);
g_K_1 = zeros(1,n_time);
g_K_2 = zeros(1,n_time);
% h
r_1 = zeros(1,n_time);  % pure
g_h_1 = zeros(1,n_time);  % uS
% syn
s_12 = zeros(1,n_time);  % pure
s_21 = zeros(1,n_time);
g_12 = zeros(1,n_time);  % uS
g_21 = zeros(1,n_time);

% set initial conditions
v_1(1) = -60;
v_2(1) = -60;
% Ca
m_Ca_1(1) = 0;
m_Ca_2(1) = 0;
g_Ca_1(1) = 0;
g_Ca_2(1) = 0;
% K
m_K_1(1) = 0;
m_K_2(1) = 0;
g_K_1(1) = 0;
g_K_2(1) = 0;
% h
r_1(1) = 0.0;
g_h_1(1) = 0;
% syn
s_12(1) = 0.0;
s_21(1) = 0.0;
g_12(1) = 0;
g_21(1) = 0;

% loop
for i = 2:n_time
    
    %%Cell 1 Ca current
    m_inf_Ca_1 = 1.0/(1.0+exp(-(v_1(i-1)-Vhalf_Ca)/Vslope_Ca));
    m_Ca_1(i) = m_inf_Ca_1+(m_Ca_1(i-1)-m_inf_Ca_1)*exp(-dt/tau_m_Ca);
    g_Ca_1(i) = GCa*m_Ca_1(i);  % uS
    %i_Ca_1 = (-g_Ca_1*(v_1(i-1)-E_Ca));
    
    %%Cell 2 Ca current
    m_inf_Ca_2 = 1.0/(1.0+exp(-(v_2(i-1)-Vhalf_Ca)/Vslope_Ca));
    m_Ca_2(i) = m_inf_Ca_2+(m_Ca_2(i-1)-m_inf_Ca_2)*exp(-dt/tau_m_Ca);
    g_Ca_2(i) = GCa*m_Ca_2(i);
    %i_Ca_2 = (-g_Ca_2*(v_2(i-1)-E_Ca));
    
    %%Cell 1 K current
    m_inf_K_1 = 1.0/(1.0+exp(-(v_1(i-1)-Vhalf_K)/Vslope_K));
    K_x_factor_1 = ((v_1(i-1)-Vhalf_K)/(2*Vslope_K));  % pure
    tau_m_K_1 = tau_0_K*2/(exp(K_x_factor_1)+exp(-(K_x_factor_1)));  % ms
    m_K_1(i) = m_inf_K_1+(m_K_1(i-1)-m_inf_K_1)*exp(-dt/tau_m_K_1);
    g_K_1(i) = GK*m_K_1(i);
    % i_K_1 = -g_K_1*(v_1(i-1)-(-80.0));
    
    %%Cell 2 K current
    m_inf_K_2 = 1.0/(1.0+exp(-(v_2(i-1)-Vhalf_K)/Vslope_K));
    K_x_factor_2 = ((v_2(i-1)-Vhalf_K)/(2*Vslope_K));
    tau_m_K_2 = tau_0_K*2/(exp(K_x_factor_2)+exp(-(K_x_factor_2)));
    m_K_2(i) = m_inf_K_2+(m_K_2(i-1)-m_inf_K_2)*exp(-dt/tau_m_K_2);
    g_K_2(i) = GK*m_K_2(i);
    %i_K_2 = -g_K_2*(v_2(i-1)-(-80.0));
    
    %%cell 1 leak current
    g_leak_1 = Gleak;  % uS
    %i_leak_1 = -g_leak_1*(v_1(i-1)-E_leak);
    
    %%cell 2 leak current
    g_leak_2 = Gleak;
    % i_leak_2 = -g_leak_2*(v_2(i-1)-E_leak);
    
    %%synapse from cell 2 to cell 1
    s_12_inf = 1/(1+exp(-(v_2(i-1)-Vhalf_syn)/Vslope_syn));
    tau_s_12 = (1-s_12_inf)*syn_tau; %% in ms
    s_12(i) = s_12_inf+(s_12(i-1)-s_12_inf)*exp(-dt/tau_s_12);
    g_12(i) = Gsyn*s_12(i); %%microS
    i_12 = -g_12(i)*(v_1(i-1)-(E_syn_12));
    
    %%synapse from 1 to 2
    s_21_inf = 1/(1+exp(-(v_1(i-1)-Vhalf_syn)/Vslope_syn));
    tau_s_21 = (1-s_21_inf)*syn_tau; % in ms
    s_21(i) = s_21_inf+(s_21(i-1)-s_21_inf)*exp(-dt/tau_s_21);
    g_21(i) = Gsyn*s_21(i); %% microS
    i_21 = -g_21(i)*(v_2(i-1)-(E_syn_21));
    
    %% Cell 1 gets h current
    %% Calculate Ih_gm1 and add to get I_12
    k_r_1 = c_r * (1+exp((v_1(i-1)-v_kr)/s_kr));
    r_inf_1 = 1/(1+exp((v_1(i-1)-v_h_half)/s_r_1));
    r_1(i) = r_inf_1 + (r_1(i-1)-r_inf_1)*exp(-dt*k_r_1);
    g_h_1(i) = Gh*r_1(i);
    %i_H_1 = g_h_1*(e_h - v_1(i-1));
    
    %% voltage update for cell 2 (no Ih)
    g_total_2 =g_K_2(i) + g_leak_2 + g_Ca_2(i) + g_21(i);
    gE_total_2=(g_K_2(i)*(E_K))+(g_leak_2*E_leak)+(g_Ca_2(i)*E_Ca)+(g_21(i)*(E_syn_21));  % nA
    v_inf_2 = (i_21 + gE_total_2)/g_total_2;
    tau_v_2 = C_m/g_total_2; %% nF/uS == ms
    v_2(i) = v_inf_2+(v_2(i-1)-v_inf_2)*exp(-dt/tau_v_2);
    
    %% voltage update for cell 1 (WITH Ih)
    g_total_1 =g_K_1(i) + g_leak_1 + g_Ca_1(i) + g_12(i) +g_h_1(i);  % uS
    gE_total_1=(g_K_1(i)*(E_K))+(g_leak_1*E_leak)+(g_Ca_1(i)*E_Ca)+(g_12(i)*(E_syn_12))+(g_h_1(i)*(e_h));  % nA
    v_inf_1 = (i_12 + gE_total_1)/g_total_1;  % mV
    tau_v_1 = C_m/g_total_1; %% nF/uS == ms
    v_1(i) = v_inf_1+(v_1(i-1)-v_inf_1)*exp(-dt/tau_v_1);
end %i-loop


if (0)
    figure (2);
    subplot(4,1,1);
    title('conductances for cell 1');
    hold on;
    plot(timevector,g_Ca_1, 'r');
    plot(timevector,g_Ca_2,':b');
    ylabel('Ca conductance (uS)');
    xlabel('time (ms)');
    legend('cell 1','cell 2');
    hold off;
    
    subplot(4,1,2);
    hold on;
    plot(timevector,g_K_1, 'r');
    plot(timevector,g_K_2,':b');
    ylabel('K conductance (uS)');
    xlabel('time (ms)');
    legend('cell 1','cell 2');
    hold off;
    
    subplot(4,1,3);
    plot(timevector,g_h_1, 'r');
    ylabel('h conductance (uS)');
    xlabel('time (ms)');
    legend('cell 1');
    
    subplot(4,1,4);
    hold on;
    plot(timevector,g_21, 'r');
    plot(timevector,g_12,':b');
    ylabel('synaptic conductance (uS)');
    xlabel('time (ms)');
    legend('cell 1','cell 2');
    hold off;
end

