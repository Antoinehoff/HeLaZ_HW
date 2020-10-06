% Run linear/nonlin simulation on a kr,kz grid and create gifs on the
% fourier and real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters and make
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
BASIC.SIMID = 'kr_kz_time_evolution'; % Name of the simulations
addpath(genpath('../matlab')) % ... add
%% PARAMETERS FOR HASEGAWA-WAKATANI SYSTEM
ALPHA = 0.1;  % Adiabadicity coefficient
KAPPA = 1.0;  % Source term (density gradient)
MU    = 0.01; % Hyper diffusivity coefficient
LK    = 1.0;  % Size of the squared frequency domain
DT    = 1e-3; % Time step
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 200;
OUTPUTS.nsave_5d = 200;
OUTPUTS.write_Ni00    = '.false.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_non_lin = '.true.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = '''results''';
OUTPUTS.rstfile0      = '''restart''';
%% Grid parameters
GRID.pmaxe = 00;  % Electron Hermite moments
GRID.jmaxe = 00;  % Electron Laguerre moments 
GRID.pmaxi = 00;  % Ion Hermite moments
GRID.jmaxi = 00;  % Ion Laguerre moments
GRID.nkr   = 128; % kr grid resolution
GRID.krmin =-LK; % kr minimum value
GRID.krmax = LK; % kr maximal value
GRID.nkz   = 128; % kz ''
GRID.kzmin =-LK; %    ''
GRID.kzmax = LK; %    ''
GRID.Pad   = 2.0; % Zero padding for dealiasing (Mx = Pad * nkx)
%% Model parameters
MODEL.CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.NON_LIN = '.true.';   % Non linear term
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
MODEL.eta_n   = -KAPPA;        % source term kappa for HW
MODEL.eta_T   = 0.0;        % Temperature
MODEL.eta_B   = 0.0;        % Magnetic
% coeff alpha for Hasegawa-Wakatani system
MODEL.lambdaD = ALPHA;
%% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
BASIC.nrun                = 100000;
BASIC.dt                  = DT;   
BASIC.tmax                = 50.0;    %time normalized to 1/omega_pe
INITIAL.RESTART           = '.true.';
INITIAL.backup_file      = '''restart''';
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile and write input file
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);
nproc = 1;
MAKE  = 'cd ..; make; cd wk';
system(MAKE);
%% Run HeLaZ
EXEC  = ' ../bin/helaz ';
RUN   = ['mpirun -np ' num2str(nproc)];
CMD   = [RUN, EXEC, INPUT];
system(CMD);
%% load results
filename = [OUTPUTS.resfile0(2:end-1),'_00.h5'];
[Nipj, p_, j_, kr, kz, time] = load_5D_data(filename, 'moments_i');
n    = squeeze(Nipj(1,1,:,:,:));
Nepj                 = load_5D_data(filename, 'moments_e');
z = squeeze(Nepj(1,1,:,:,:));
phi                  = load_2D_data(filename, 'phi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
if strcmp(MODEL.NON_LIN,'.true.'); LINEARITY = 'nl';
else; LINEARITY = 'lin'; end;
default_plots_options
SAVEFIG = 1;
default_plots_options
% Variables
xmax = GRID.nkr/(2*(GRID.krmax-GRID.krmin));
xmin = -xmax;
ymax = GRID.nkz/(2*(GRID.kzmax-GRID.kzmin));
ymin = -ymax;
x = linspace(xmin,xmax,GRID.nkr);
y = linspace(ymin,ymax,GRID.nkz);
disp('');
skip    = 1;    % To skip some frames
delay   = 0.03; % Speed of the movie (smaller for faster)
Fn      = zeros(size(n));
Fz      = zeros(size(n));
Fp      = zeros(size(n));
for it = 1:numel(time)
    F_         = real(ifft2(n(:,:,it),floor(GRID.nkr),floor(GRID.nkz)));
    Fn(:,:,it) = ifftshift(F_(1:GRID.nkr, 1:GRID.nkr));
    F_         = real(ifft2(z(:,:,it),floor(GRID.nkr),floor(GRID.nkz)));
    Fz(:,:,it) = ifftshift(F_(1:GRID.nkr, 1:GRID.nkr));
    F_         = real(ifft2(phi(:,:,it),floor(GRID.nkr),floor(GRID.nkz)));
    Fp(:,:,it) = ifftshift(F_(1:GRID.nkr, 1:GRID.nkr));
end
%% Show last frame of n, zeta and phi
it = numel(time);
fig = figure;
    subplot(231)
        pclr = pcolor(kr,kz,log(abs(n(:,:,it))));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('kr');ylabel('kz');
        title(['log$|n|$, $t \approx$', sprintf('%.3d',ceil(time(it)))])
    subplot(234)
        pclr = pcolor(x,y,Fn(:,:,end));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('x');ylabel('y');
        title(['$\mathcal{F}^{-1}\{n\}$, $t \approx$', sprintf('%.3d',ceil(time(it)))])
    subplot(232)
        pclr = pcolor(kr,kz,log(abs(z(:,:,it))));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('kr');ylabel('kz');
        title(['log$|\zeta|$, $t \approx$', sprintf('%.3d',ceil(time(it)))])
    subplot(235)
        pclr = pcolor(x,y,Fz(:,:,end));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('x');ylabel('y');
        title(['$\mathcal{F}^{-1}\{\zeta\}$, $t \approx$', sprintf('%.3d',ceil(time(it)))])
    subplot(233)
        pclr = pcolor(kr,kz,log(abs(phi(:,:,it))));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('kr');ylabel('kz');
        title(['log$|\phi|$, $t \approx$', sprintf('%.3d',ceil(time(it)))])
    subplot(236)
        pclr = pcolor(x,y,Fp(:,:,end));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('x');ylabel('y');
        title(['$\mathcal{F}^{-1}\{\phi\}$, $t \approx$', sprintf('%.3d',ceil(time(it)))])
FIGNAME = [LINEARITY,'-',num2str(GRID.nkr),'x',num2str(GRID.nkz),...
           '-dt-',num2str(BASIC.dt),'_last_frame'];
save_figure
%% Gif creation
%% moments in frequency space
GIFNAME = [LINEARITY,'-',num2str(GRID.nkr),'x',num2str(GRID.nkz),...
           '-dt-',num2str(BASIC.dt),'-Ni00'];
X = kr; Y = kz; T = time; FIELD = log(abs(n)); FIELDNAME = 'log$|n|$';
SATURATE = 0; DELAY = 0.1;
%create_gif
%% moments in real space
GIFNAME = [LINEARITY,'-',num2str(GRID.nkr),'x',num2str(GRID.nkz),...
           '-dt-',num2str(BASIC.dt),'-Mi00'];
X = kr; Y = kz; T = time; FIELD = Fn; FIELDNAME = '$\mathcal{F}^{-1}\{N_i^{00}\}$';
SATURATE = 0; DELAY = 0.1;
%create_gif
%% non linear term
if strcmp(OUTPUTS.write_non_lin,'.true.')
Sipj = load_5D_data(filename, 'Sipj');
Si00    = squeeze(Sipj(1,1,:,:,1:skip:end));
GIFNAME = [LINEARITY,'-',num2str(GRID.nkr),'x',num2str(GRID.nkz),...
           '-dt-',num2str(BASIC.dt),'-Si00'];
X = kr; Y = kz; T = time; FIELD = log(abs(Si00)); FIELDNAME = 'log$|S_i^{00}|$';
SATURATE = 0; DELAY = 0.1;
%create_gif
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%