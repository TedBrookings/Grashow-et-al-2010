function Currents = SimDynamicClampCurrents(VBio, VModel, dt)

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
E_syn_modbio = -90; % mV, from mod to bio, should  be -90
E_syn_21 = -70; % mV, from bio to mod should be -70
GCa = 200e-3;  % uS
GK = 200e-3;  % uS

% Model Ca current
m_inf_Ca = 1.0/(1.0+exp(-(V_mod-Vhalf_Ca)/Vslope_Ca));
m_Ca = m_inf_Ca + (m_Ca - m_inf_Ca) * exp(-dt_sim / tau_m_Ca);
g_Ca_mod = gmax_Ca * m_Ca;
i_Ca = (-g_Ca_mod * (V_mod - E_Ca));

%Model K current
m_inf_K = 1.0/(1.0+exp(-(V_mod-Vhalf_K)/Vslope_K));
K_x_factor = ((V_mod - Vhalf_K[k]) / (2 * Vslope_K));
tau_m_K = tau_0_K*2/(exp(K_x_factor)+exp(-(K_x_factor)));
m_K = m_inf_K+(m_K-m_inf_K)*exp(-dt_sim/tau_m_K);
g_K_mod = gmax_K*m_K;
i_K_mod = -g_K_mod*(V_mod-(-80.0));

%Model leak current
g_leak_mod = gmax_leak;
i_leak_mod = -g_leak_mod*(V_mod-E_leak);

%synapse from biological to model (from 1 to 2, syn_21)
s_modbio_inf = 1.0 / (1.0 + exp(-(VBio - Vhalf_syn) / Vslope_syn));
tau_s_modbio = (1.0 - s_modbio_inf) * syn_tau; % in ms
s_modbio = s_modbio_inf + (s_modbio - s_modbio_inf) * exp(-dt / tau_s_modbio);
g_modbio = g_syn * s_modbio; %microS
i_modbio = -g_modbio * (VModel - (-90.0));


%synapse from model to biological (from 2 to 1, syn_12)
s_biomod_inf = 1.0 / (1.0 + exp(-(VModel - Vhalf_syn) / Vslope_syn));
tau_s_biomod = (1.0 - s_biomod_inf) * syn_tau; % in ms
s_biomod = s_biomod_inf + (s_biomod - s_biomod_inf) * exp(-dt / tau_s_biomod);
g_biomod = g_syn * s_biomod; %microS
i_biomod = -g_biomod * (VBio - (-90.0));

  
% Calculate Ih_gm1 and add to get I_12
k_r_1 = c_r * (1.0 + exp((VBio - v_kr) / s_kr));
r_inf_1 = 1.0 / (1.0 + exp((Vbio - v_h_half) / s_r_1));
r_bio = r_inf_1 + (r_bio - r_inf_1) * exp(-dt * k_r_1);
i_H_bio = (g_h * r_bio) * (e_h - VBio);
i_to_bio = i_biomod + i_H_bio;

% voltage update
g_total_mod = g_K_mod + g_leak_mod + g_Ca_mod + g_modbio;
gE_total = (g_K_mod*(-80)) + (g_leak_mod * E_leak) ...
    + (g_Ca_mod * E_Ca) + (g_modbio * (-70.0));  % nA
v_inf_mod = (mod_DC_I + i_modbio + gE_total) / g_total_mod;
v_tau_mod = C_m / g_total_mod;  % nF/uS == ms
V_mod = v_inf_mod + (V_mod - v_inf_mod) * exp(-dt_sim / v_tau_mod);

return
