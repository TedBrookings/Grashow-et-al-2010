function SimData = SimDynamicClamp(t, VBio, VModel, g_syn, g_h, SynapseForm)
if(nargin < 6)
  SynapseForm = 0;
end
%VBio = round(2 * VBio) / 2;

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
gmax_leak =200e-3; % uS
E_leak = -60;  % mV
Vhalf_syn = -45;  % mV
Vslope_syn = 5;  % mV
E_K = -80;  % mV
E_syn_modbio = -90; % mV, from mod to bio, should  be -90
E_syn_21 = -70; % mV, from bio to mod should be -70
gmax_Ca = 200e-3;  % uS
gmax_K = 200e-3;  % uS

gmax_H_bio = 0.001 * g_h;
gmax_biomod = 0.001 * g_syn;
gmax_modbio = 0.001 * g_syn;
Vhalf_syn_modbio = Vhalf_syn;
Vhalf_syn_biomod = Vhalf_syn;
Vslope_syn_modbio = Vslope_syn;
Vslope_syn_biomod = Vslope_syn;

mod_DC_I = 0;

dt = 200e-3;  %200 us, expressed in ms
fprintf('dt = %g\n', dt)
dt_sim = dt/10.0;
t_n = 2 * t(1) - t(2);

%initializing state variables
m_Ca = 0.0;
m_K = 0.0;
V_mod = -60.0;
m_H_bio = 0.0;
s_biomod = 0.0;
s_modbio = 0.0;
r_bio = 0.0;

%initialize dV_mod
NumT = length(t);
V_mod_calc = zeros(size(VModel));
I_modbio = zeros(size(VModel)); 

V_mod_calc(1) = V_mod;
for n = 1:NumT
  % for each timestep
  %V_mod = VModel(n);
  V_bio = VBio(n);
  V_mod0 = V_mod;
  
  while(t_n < t(n))
    for i=0:9
      %Model Ca current
      m_inf_Ca = 1.0 / (1.0 + exp(-(V_mod - Vhalf_Ca) / Vslope_Ca));
      m_Ca = m_inf_Ca + (m_Ca - m_inf_Ca) * exp(-dt_sim / tau_m_Ca);
      g_Ca_mod = gmax_Ca * m_Ca;
      i_Ca = (-g_Ca_mod * (V_mod - E_Ca));
      
      %Model K current
      m_inf_K = 1.0 / (1.0 + exp(-(V_mod - Vhalf_K) / Vslope_K));
      K_x_factor = ((V_mod - Vhalf_K) / (2 * Vslope_K));
      tau_m_K = tau_0_K * 2 / (exp(K_x_factor) + exp(-(K_x_factor)));
      m_K = m_inf_K + (m_K - m_inf_K) * exp(-dt_sim / tau_m_K);
      g_K_mod = gmax_K * m_K;
      i_K_mod = -g_K_mod * (V_mod - (-80.0));
      
      %Model leak current
      g_leak_mod = gmax_leak;
      i_leak_mod = -g_leak_mod * (V_mod - E_leak);
      
      %synapse from biological to model (from 1 to 2, syn_21)
      s_modbio_inf = 1.0 / (1.0 + exp(-(V_bio - Vhalf_syn_modbio) ...
				      / Vslope_syn_modbio));
      tau_s_modbio = (1 - s_modbio_inf) * syn_tau; % in ms
      s_modbio = s_modbio_inf + (s_modbio - s_modbio_inf) ...
                 * exp(-dt_sim / tau_s_modbio);
      g_modbio = gmax_modbio * s_modbio; %microS
      switch(SynapseForm)
        case 0,  %as it was
          i_modbio = -g_modbio * (V_mod - (-90.0));
      
          % voltage update
          g_total_mod =g_K_mod + g_leak_mod + g_Ca_mod + g_modbio;
          gE_total = (g_K_mod * (-80)) + (g_leak_mod * E_leak) ...
	  + (g_Ca_mod * E_Ca) + (g_modbio * (-70.0));  % nA
          v_inf_mod = (mod_DC_I + i_modbio + gE_total) / g_total_mod;
        case 1, %with i_modbio removed
          i_modbio = -g_modbio * (V_mod - (-70.0));
      
          % voltage update
          g_total_mod =g_K_mod + g_leak_mod + g_Ca_mod + g_modbio;
          gE_total = (g_K_mod * (-80)) + (g_leak_mod * E_leak) ...
	  + (g_Ca_mod * E_Ca) + (g_modbio * (-70.0));  % nA
          v_inf_mod = (mod_DC_I + gE_total) / g_total_mod;        
        case 2, %with i_modbio added according to Abbot scheme
          i_modbio = -g_modbio * (2 * V_mod - (-90.0 - 70.0));
      
          % voltage update
          g_total_mod =g_K_mod + g_leak_mod + g_Ca_mod + 2* g_modbio;
          gE_total = (g_K_mod * (-80)) + (g_leak_mod * E_leak) ...
	  + (g_Ca_mod * E_Ca) + (g_modbio * (-70.0 - 90.0));  % nA
          v_inf_mod = (mod_DC_I + gE_total) / g_total_mod;        
      end
      v_tau_mod = C_m / g_total_mod;  % nF/uS == ms
      V_mod = v_inf_mod + (V_mod - v_inf_mod) * exp(-dt_sim / v_tau_mod);
    end
    V_mod = VModel(n);
    
    % synapse from model to biological (from 2 to 1, syn_12)
    s_biomod_inf = 1 / (1 + exp(-(V_mod - Vhalf_syn_biomod) ...
				/ Vslope_syn_biomod));
    tau_s_biomod = (1 - s_biomod_inf) * syn_tau; % in ms
    s_biomod = s_biomod_inf ...
	+ (s_biomod - s_biomod_inf) * exp(-dt/tau_s_biomod);
    g_biomod = gmax_biomod * s_biomod; %microS
    i_biomod = -g_biomod * (V_bio - (-90));
    
    % Other h current
    % Calculate Ih_gm1 and add to get I_12
    k_r_1 = c_r * (1+exp((V_bio-v_kr)/s_kr));
    r_inf_1 = 1/(1+exp((V_bio-v_h_half)/s_r_1));
    r_bio = r_inf_1 + (r_bio-r_inf_1)*exp(-dt*k_r_1);
    i_H_bio = (gmax_H_bio *r_bio)*(e_h - V_bio);
    i_to_bio = i_biomod+i_H_bio;
    
    t_n = t_n + dt;
    %V_mod = round(2 * V_mod) / 2;
  end
  
  %My little addition:  record the change in model voltage
  %dV_mod(n) = V_mod - V_mod0;
  V_mod_calc(n) = V_mod;
  I_modbio(n) = i_to_bio; 
end

SimData.VModelCalc = V_mod_calc;
SimData.I_InjectedCalc = I_modbio;
return


