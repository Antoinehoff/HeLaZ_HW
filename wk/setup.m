clear all;
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS FOR HASEGAWA-WAKATANI SYSTEM
ALPHA   = 1.0;    % Adiabadicity coefficient
KAPPA   = 1.0;    % Source term (density gradient)
MU      = 2e-4;   % Hyper diffusivity coefficient
%% GRID PARAMETERS
N       = 64;    % Frequency gridpoints
L       = 40;     % Size of the squared frequency domain
%% TIME PARAMETERS 
TMAX    = 100;  % Maximal time unit
DT      = 1e-3;   % Time step
SPS     = 10;      % Sampling per time unit
RESTART = 0;     % To restart from last checkpoint
%% OPTIONS
SIMID   = '';  % Name of the simulation
HALF    = 1; % To solve only on half a domain
NON_LIN = 1; % activate non-linearity 
MHW     = 0; % Cancel non zonal term in the adiabadic response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ________________________________________________________________________
% Naming and creating input file
params   = ['a_',num2str(ALPHA),'_mu_%0.0e_'];
params   = sprintf(params,MU);
if ~NON_LIN; params = ['lin_',params]; end;
if MHW;    model  = '_MHW_'; else; model = '_HW_'; end;
resolution = [num2str(N),'x',num2str(N*(1-HALF/2)),'_'];
gridname   = ['L_',num2str(L),'_'];
BASIC.SIMID = [SIMID,model,resolution,gridname,params];
BASIC.nrun       = 1e8;
BASIC.dt         = DT;   
BASIC.tmax       = TMAX;    %time normalized to 1/omega_pe
if RESTART; BASIC.RESTART = '.true.'; else; BASIC.RESTART = '.false.';end;
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = floor(1.0/SPS/DT);
OUTPUTS.nsave_5d = floor(1.0/SPS/DT);
OUTPUTS.write_Ni00    = '.false.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_non_lin = '.true.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = ['''',BASIC.SIMID,''''];
OUTPUTS.rstfile0      = ['''','../checkpoint/cp_',BASIC.SIMID,''''];
% Grid parameters
GRID.pmaxe = 00;  % Electron Hermite moments
GRID.jmaxe = 00;  % Electron Laguerre moments 
GRID.pmaxi = 00;  % Ion Hermite moments
GRID.jmaxi = 00;  % Ion Laguerre moments
GRID.Nr    = N; % r grid resolution
GRID.Lr    = L; % r length
GRID.Nz    = N; % z ''
GRID.Lz    = L; % z ''
% Model parameters
MODEL.CO      = 0;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
if NON_LIN; MODEL.NON_LIN = '.true.'; else; MODEL.NON_LIN = '.false.';end;
if MHW; MODEL.MHW = '.true.'; else; MODEL.MHW = '.false.';end;
MODEL.nu      = MU; % hyper diffusive coefficient nu for HW
% temperature ratio T_a/T_e
MODEL.tau_e   = 1.0;
MODEL.tau_i   = 1.0;
% mass ratio sqrt(m_a/m_i)
MODEL.sigma_e = 0.0233380;
MODEL.sigma_i = 1.0;
% charge q_a/e
MODEL.q_e     =-1.0;
MODEL.q_i     = 1.0;
% gradients L_perp/L_x
MODEL.eta_n   = KAPPA;        % source term kappa for HW
MODEL.eta_T   = 0.0;        % Temperature
MODEL.eta_B   = 0.0;        % Magnetic
% coeff alpha for Hasegawa-Wakatani system
MODEL.lambdaD = ALPHA;
% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
INITIAL.only_Na00         = '.false.';
INITIAL.initback_moments  = 1.0e-4;
INITIAL.initnoise_moments = 5.0e-5;
INITIAL.iseed             = 42;
INITIAL.selfmat_file = ...
    ['''../iCa/self_Coll_GKE_0_GKI_0_ESELF_1_ISELF_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10.h5'''];
INITIAL.eimat_file = ...
    ['''../iCa/ei_Coll_GKE_0_GKI_0_ETEST_1_EBACK_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];
INITIAL.iemat_file = ...
    ['''../iCa/ie_Coll_GKE_0_GKI_0_ITEST_1_IBACK_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];

%% Compile and write input file
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);
nproc = 1;
MAKE  = 'cd ..; make; cd wk';
system(MAKE);
%%
disp(['Set up ', BASIC.SIMID]);