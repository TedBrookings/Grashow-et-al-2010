// Program to connect a Morris Lecar neuron to an GM neuron- modified from LP_ML_hybrid program. Rachel G, 3/2/09

#include <model.h>
#include <math.h>

// Shared memory direct values-entered by user
static double gmax_Ca[MAX_FUNCTIONS];
static double Vhalf_Ca[MAX_FUNCTIONS];
static double Vslope_Ca[MAX_FUNCTIONS];
static double gmax_K[MAX_FUNCTIONS];
static double Vhalf_K[MAX_FUNCTIONS];
static double Vslope_K[MAX_FUNCTIONS];
static double gmax_H_bio[MAX_FUNCTIONS];
static double gmax_leak[MAX_FUNCTIONS];
static double E_leak[MAX_FUNCTIONS];
static double gain[MAX_FUNCTIONS];
static double Vhalf_syn_modbio[MAX_FUNCTIONS];
static double Vslope_syn_modbio[MAX_FUNCTIONS];
static double gmax_biomod[MAX_FUNCTIONS];
static double gmax_modbio[MAX_FUNCTIONS];
static double mod_DC_I[MAX_FUNCTIONS];
static double Vhalf_syn_biomod[MAX_FUNCTIONS];
static double Vslope_syn_biomod[MAX_FUNCTIONS];

// Shared memory states- updated at every time step
static double V_mod[MAX_FUNCTIONS];
static double m_Ca[MAX_FUNCTIONS];
static double m_K[MAX_FUNCTIONS];
static double m_H_bio[MAX_FUNCTIONS];
static double s_biomod[MAX_FUNCTIONS];
static double s_modbio[MAX_FUNCTIONS];
static double r_bio[MAX_FUNCTIONS];

////Added by RG to make program new RTLDC compatible-3/2/09
extern double pos(double x)
    {
    return (x>=0)?x:0.0;
    }

extern void GM_ML_hybrid_main (double *data,int mem_index)
  {
      
  static int i,k;
  static double dt,gE_total,i_Ca,V_bio,m_inf_Ca; 
  static double g_Ca_mod,m_inf_K,tau_m_K, K_x_factor;
  static double g_K_mod,i_K_mod,g_total_mod,tau_s_modbio;
  static double s_biomod_inf,s_modbio_inf,v_tau_mod,tau_s_biomod;
  static double g_leak_mod,i_leak_mod,v_inf_mod, i_to_bio;
  static double dt_sim, g_biomod, i_biomod,  g_modbio, i_modbio;
  static double C_m, E_Ca, tau_m_Ca, tau_0_K, syn_tau;
  static double s_r_1,v_h_half,c_r, v_kr, e_h, s_kr, k_r_1;
  static double r_inf_1, i_H_bio, g_H_bio; 

  k = mem_index;
  dt = shmem->dT_ms; //ms
  dt_sim = dt/10.0;  // ms
  V_bio = input(0)*1e3; //mV
  C_m = 10.0;  // cell capacitance, nF
  E_Ca = +100.0; //mV
  tau_m_Ca = 0.001;  // ms
  tau_0_K = 1.0/0.002;  // ms
  syn_tau = 100; //ms
  s_r_1 = 7;
  v_h_half = -50;
  c_r= .00033; //1/ms
  v_kr = -110; //mV
  e_h = -10; //mV
  s_kr = -13; //mV



  // for each timestep
  for(i=0;i<10;i++)
    {
    //Model Ca current
    m_inf_Ca = 1.0/(1.0+exp(-(V_mod[k]-Vhalf_Ca[k])/Vslope_Ca[k]));
    m_Ca[k] = m_inf_Ca+(m_Ca[k]-m_inf_Ca)*exp(-dt_sim/tau_m_Ca);
    g_Ca_mod = gmax_Ca[k]*m_Ca[k];
    i_Ca = (-g_Ca_mod*(V_mod[k]-E_Ca));

    //Model K current
    m_inf_K = 1.0/(1.0+exp(-(V_mod[k]-Vhalf_K[k])/Vslope_K[k]));
    K_x_factor = ((V_mod[k]-Vhalf_K[k])/(2*Vslope_K[k]));
    tau_m_K = tau_0_K*2/(exp(K_x_factor)+exp(-(K_x_factor)));
    m_K[k] = m_inf_K+(m_K[k]-m_inf_K)*exp(-dt_sim/tau_m_K);
    g_K_mod = gmax_K[k]*m_K[k];
    i_K_mod = -g_K_mod*(V_mod[k]-(-80.0));

    //Model leak current
    g_leak_mod = gmax_leak[k];
    i_leak_mod = -g_leak_mod*(V_mod[k]-E_leak[k]);
	
    //synapse from biological to model (from 1 to 2, syn_21)
    s_modbio_inf = 1/(1+exp(-(V_bio-Vhalf_syn_modbio[k])/Vslope_syn_modbio[k]));
    tau_s_modbio = (1-s_modbio_inf)*syn_tau; // in ms
    s_modbio[k] = s_modbio_inf+(s_modbio[k]-s_modbio_inf)*exp(-dt_sim/tau_s_modbio);
    g_modbio = gmax_modbio[k]*s_modbio[k]; //microS
    i_modbio = -g_modbio*(V_mod[k]-(-90.0));
  
    // voltage update
    g_total_mod =g_K_mod + g_leak_mod + g_Ca_mod + g_modbio;
    gE_total=(g_K_mod*(-80))+(g_leak_mod*E_leak[k])+(g_Ca_mod*E_Ca)+(g_modbio*(-70.0));  // nA
    v_inf_mod = (mod_DC_I[k] + i_modbio + gE_total)/g_total_mod;
    v_tau_mod = C_m/g_total_mod;  // nF/uS == ms
    V_mod[k] = v_inf_mod+(V_mod[k]-v_inf_mod)*exp(-dt_sim/v_tau_mod);
    }
 
  //synapse from model to biological (from 2 to 1, syn_12)
  s_biomod_inf = 1/(1+exp(-(V_mod[k]-Vhalf_syn_biomod[k])/Vslope_syn_biomod[k]));
  tau_s_biomod = (1-s_biomod_inf)*syn_tau;// in ms
  s_biomod[k] = s_biomod_inf+(s_biomod[k]-s_biomod_inf)*exp(-dt/tau_s_biomod);
  g_biomod = gmax_biomod[k]*s_biomod[k]; //microS
  i_biomod = -g_biomod*(V_bio-(-90));

  //Other h current
  // Calculate Ih_gm1 and add to get I_12
  k_r_1 = c_r * (1+exp((V_bio-v_kr)/s_kr));
  r_inf_1 = 1/(1+exp((V_bio-v_h_half)/s_r_1));
  r_bio[k] = r_inf_1 + (r_bio[k]-r_inf_1)*exp(-dt*k_r_1);
  i_H_bio = (gmax_H_bio[k] *r_bio[k])*(e_h - V_bio);
  i_to_bio = i_biomod+i_H_bio;


  //can watch during experiment with drop down menu
  state(0) = V_bio; //mV
  state(1) = V_mod[k]; //mV
  state(2) = i_to_bio; //nA
  state(3) = i_modbio; //nA
  state(4) = g_biomod*1000.0; //uS->nS
  state(5) = g_modbio*1000.0; //uS->nS
  state(6) = g_H_bio*1000.0; //uS->nS

  output(0) = output(0) + i_to_bio/gain[k]; //V
  output(1) = output(1) + (V_mod[k]*0.01); //V
  }

extern void GM_ML_hybrid_update(double *data,int mem_index,int first_time)
  {
  static int k;
  k = mem_index;

  //initial values for direct shared memory states
  if (first_time) 
    {
    shmem->parameter[k][ 0] = 200.0;  // dflt gmax_Ca, nS
    shmem->parameter[k][ 1] =   -20.0;  // dflt Vhalf_Ca, mV
    shmem->parameter[k][ 2] =   15;  // dflt Vslope_ca, mV
    shmem->parameter[k][ 3] = 200.0;  // dflt gmax_K, nS 
    shmem->parameter[k][ 4] =  -20.0;  // dflt Vhalf_K, mV
    shmem->parameter[k][ 5] =   15;  // dflt Vslope_K,mV
    shmem->parameter[k][ 6] =   0.0;  // dflt biological H conductance
    shmem->parameter[k][ 7] =  200.0;  // dflt leak conductance, nS
    shmem->parameter[k][ 8] = -60.0; // dflt leak reversal potential, mV
    shmem->parameter[k][9] =   1.0; // dflt headstage gain
    shmem->parameter[k][10] = -25.0;  // dflt synapse half activation for biological, mV
    shmem->parameter[k][11] =   15.0;  // dflt synapse slope factor for biological, mV
    shmem->parameter[k][12] =   0.0;  // dflt mod to bio synapse, nS
    shmem->parameter[k][13] =   0.0;  // dflt bio to mod synapse, nS
    shmem->parameter[k][14] =   0.0; //dflt for DC current into model,nA
    shmem->parameter[k][15] =   -25; //dflt synapse half activation for model, mV
    shmem->parameter[k][16] =   15; //dflt synapse slope factor for model, mV
    }

  //getting variables from shared memory and putting into variables that I will use for my program
  //should be the same number as the list above
  gmax_Ca[k] = shmem->parameter[k][0]/1000.0; //gmax_ca, nS->uS
  Vhalf_Ca[k] = shmem->parameter[k][1]; //Vhalf_Ca, mV
  Vslope_Ca[k] = shmem->parameter[k][2]; //Vslope_Ca, mV
  gmax_K[k] = shmem->parameter[k][3]/1000.0; //gmax_K, nS->uS
  Vhalf_K[k] = shmem->parameter[k][4]; //Vhalf_K, mV
  Vslope_K[k] = shmem->parameter[k][5]; //Vslope_K, mV
  gmax_H_bio[k] = shmem->parameter[k][6]/1000.0;//gmax_H_bio, nS->uS
  gmax_leak[k] = shmem->parameter[k][7]/1000.0; //gmax_leak, nS->uS
  E_leak[k] = shmem->parameter[k][8]; // reversal potential for leak conductance
  gain[k] = shmem->parameter[k][9]*10.0; //headstage gain *10, current to voltage command scaling factor
  Vhalf_syn_modbio[k] = shmem->parameter[k][10]; // half activation synaptic voltage for biological, mV
  Vslope_syn_modbio[k] = shmem->parameter[k][11]; //synaptic slope factor for biological, mV
  gmax_biomod[k] = shmem->parameter[k][12]/1000.0; // model to bio maximal conductance, nS->uS
  gmax_modbio[k] = shmem->parameter[k][13]/1000.0; // model to bio maximal condunctance, nS->uS
  mod_DC_I[k] = shmem->parameter[k][14]; // DC current into model neuron
  Vhalf_syn_biomod[k] = shmem->parameter[k][15];// synapse half activation for model, mV
  Vslope_syn_biomod[k] = shmem->parameter[k][16]; // synapse slope factor for model, mV 

  // Reset state variables- for every gating variable
  //if I change parameters, I will want to reset the gating variables
  V_mod[k] = -60.0;  // mV
  m_Ca[k] =  0.0; 
  m_K[k] =  0.0;
  m_H_bio[k] =  0.0; 
  s_biomod[k] = 0.0;
  s_modbio[k] = 0.0;

  }

int init_module(void)
  {

#if defined RTLINUX
  shmem = (rtldc_shm *)mbuff_alloc(SHM_ADDR,sizeof(rtldc_shm));
#elif defined RTAI
  shmem = (rtldc_shm *)rtai_kmalloc(nam2num(SHM_ADDR),sizeof(rtldc_shm));
#else
#error "Must #define RTLINUX or RTAI"
#endif

  shmem->main_fxn    = &GM_ML_hybrid_main;
  shmem->update_fxn  = &GM_ML_hybrid_update;
  
  return 0;
  }

void cleanup_module(void)
  {

#if defined RTLINUX
  mbuff_free (SHM_ADDR,(void *)shmem); 
#elif defined RTAI
  //rtai_kfree((int)shmem);
  rtai_kfree(nam2num(SHM_ADDR));
#else
#error "Must #define RTLINUX or RTAI"
#endif

  }
