%% load results
SID = 00;
filename = [BASIC.SIMID,'_','%.2d.h5'];
filename = sprintf(filename,SID); disp(['Analysing ',filename])
[Nipj, p_, j_, kr, kz, Ts] = load_5D_data(filename, 'moments_i');
Nepj                       = load_5D_data(filename, 'moments_e');
NN      = squeeze(Nipj(1,1,:,:,:));
ZZ      = squeeze(Nepj(1,1,:,:,:));
PP      = load_2D_data(filename, 'phi');
Ts      = Ts';
Ns      = numel(Ts);
dt      = mean(diff(Ts));
if strcmp(OUTPUTS.write_non_lin,'.true.')
    Sipj    = load_5D_data(filename, 'Sipj');
    Sepj    = load_5D_data(filename, 'Sepj');
    SN      = squeeze(Sipj(1,1,:,:,:));
    SZ      = squeeze(Sepj(1,1,:,:,:));
end
%% Build grids
Nkr = numel(kr); Nkz = numel(kz);
[KZ,KR] = meshgrid(kz,kr);
Lkr = max(kr)-min(kr); Lkz = max(kz)-min(kz);
dkr = Lkr/(Nkr-1); dkz = Lkz/(Nkz-1);
Lk = max(Lkr,Lkz);

dr = 2*pi/Lk; dz = 2*pi/Lk;
Nr = max(Nkr,Nkz);         Nz = Nr;
r = dr*(-Nr/2:(Nr/2-1)); Lr = max(r)-min(r);
z = dz*(-Nz/2:(Nz/2-1)); Lz = max(z)-min(z);
[YY,XX] = meshgrid(z,r);
%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IFFT
zeta   = zeros(Nr,Nz);
ni     = zeros(Nr,Nz);
phi    = zeros(Nr,Nz);

for it = 1:numel(PP(1,1,:))
    ZZ_ = ZZ(:,:,it); NN_ = NN(:,:,it); PP_ = PP(:,:,it);
    F_          = (ifft2((ZZ_),Nr,Nz));
    zeta(:,:,it)= real(fftshift(F_));
    F_          = (ifft2((NN_),Nr,Nz));
    ni(:,:,it)  = real(fftshift(F_));
    F_          = (ifft2((PP_),Nr,Nz));
    phi(:,:,it) = real(fftshift(F_));
end

%% Post processing
phi_ST = zeros(Nz,Ns);   % Space-Time diagram of ES potential
n_ST   = zeros(Nz,Ns);   % Space-Time diagram of densitz
phi_00 = zeros(1,Ns);    % Time evolution of ES potential at origin
ni_00  = zeros(1,Ns);    % Time evolution of density at origin
NN_11  = zeros(1,Ns);    % Time evolution of F density at 1,1
Sn_norm= zeros(1,Ns);    % Time evolution of the amp of density nonlin term
Sz_norm= zeros(1,Ns);    % Time evolution of the amp of vorti. nonlin term
E_pot  = zeros(1,Ns);    % Potential energy n^2
E_kin  = zeros(1,Ns);    % Kinetic energy grad(phi)^2
W      = zeros(1,Ns);    % Enstrophy
G_n    = zeros(1,Ns);    % Background density energy source
G_a    = zeros(1,Ns);    % Adiabadicity energy sink
D_E    = zeros(1,Ns);    % Dissipative energy term
D_W    = zeros(1,Ns);    % Dissipative vorticity term
ExB    = zeros(1,Ns);    % ExB drift intensity \propto |\grad \phi|
CFL    = zeros(1,Ns);    % CFL time step
Ddr = 1i*KR; Ddz = 1i*KZ; lapl   = Ddr.^2 + Ddz.^2; 
[~, ikr1] = min(abs(kr-1));
[~, ikz1] = min(abs(kz-1));
for it = 1:numel(PP(1,1,:))
    ZZ_ = ZZ(:,:,it); NN_ = NN(:,:,it); PP_ = PP(:,:,it);
    phi_ST(:,it)= phi(:,z==0,it);
    n_ST  (:,it)= ni(:,z==0,it);
    phi_00(it)  = phi(r==0,z==0,it);
    ni_00(it)   = ni(r==0,z==0,it);
    NN_11(it)   = NN(ikr1,ikz1,it);
    Sn_norm(it) = sum(sum(abs(SN(:,:,it))));
    Sz_norm(it) = sum(sum(abs(SZ(:,:,it))));
    E_pot(it)   = pi/Lr/Lz*sum(sum(abs(NN_).^2))/Nkr/Nkr; % integrate through Parseval id
    E_kin(it)   = pi/Lr/Lz*sum(sum(abs(Ddr.*PP_).^2+abs(Ddz.*PP_).^2))/Nkr/Nkr;
    W(it)       = pi/Lr/Lz*sum(sum(abs(NN_ - lapl.*PP_).^2))/Nkr/Nkr;
    G_n(it)     =real(-2*pi*KAPPA/Lr/Lz*sum(sum((NN_.*conj(Ddz.*PP_))))/Nkr/Nkr);
    G_a(it)     = 2*pi*ALPHA/Lr/Lz*sum(sum(abs(NN_-PP_).^2))/Nkr/Nkr;
    D_E(it)     = 2*pi*MU   /Lr/Lz*sum(sum(abs(lapl.*NN_).^2 + abs(Ddr.*(lapl.*PP_)).^2 + abs(Ddz.*(lapl.*PP_)).^2))/Nkr/Nkr;
    D_W(it)     = real(2*pi*MU   /Lr/Lz*sum(sum(abs(lapl.*NN_).^2 + abs(lapl.^2.*PP_).^2 - 2*(lapl.*NN_).*conj(lapl.^2.*PP_)))/Nkr/Nkr);
    ExB(it)     = max(max(max(abs(phi(3:end,:,it)-phi(1:end-2,:,it))/(2*dr))),max(max(abs(phi(:,3:end,it)-phi(:,1:end-2,it))'/(2*dz))));
    CFL(it)     = 4*min([dr^2/MU,dr/KAPPA,ExB(it)]);
end
E_kin_KZ = mean(mean(abs(Ddr.*PP(:,:,it)).^2+abs(Ddz.*PP(:,:,it)).^2,3),1);
E_kin_KR = mean(mean(abs(Ddr.*PP(:,:,it)).^2+abs(Ddz.*PP(:,:,it)).^2,3),2);
dEdt     = diff(E_pot+E_kin)./diff(Ts);
dWdt     = diff(W)./diff(Ts);
%% PLOTS
%% Time evolutions
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',SID)];
    subplot(221)
        semilogy(Ts,abs(ni_00),'-','DisplayName','$n$')
        grid on; xlabel('$t$'); ylabel('$|n(x=0,y=0)|$');
    subplot(222)
        semilogy(Ts,abs(NN_11),'-','DisplayName','$\phi$')
        grid on; xlabel('$t$'); ylabel('$|\tilde n(k_r\approx 1,k_z\approx 1)|$');
    subplot(223)
        semilogy(Ts,E_kin+E_pot,'-','DisplayName','$\sum|ik\tilde\phi_i|^2+\sum|\tilde n_i|^2$')
        hold on;
        if ~NON_LIN % Plot linear growth rate
            [Gmax, GG] = HW_lin_disp_rel(ALPHA,KAPPA,MU,KR,KZ);
            semilogy(Ts,(E_pot(end)+E_kin(end)).*exp(2*Gmax.*(Ts-Ts(end))),'--k','DisplayName','$\exp(2\gamma_{\max}t)$')
        end
        grid on; xlabel('$t$'); ylabel('$E$'); legend('show');
    subplot(224)
        semilogy(Ts,Sn_norm,'-','DisplayName','$\sum|S_n|$'); 
        hold on;
        semilogy(Ts,Sz_norm,'-','DisplayName','$\sum|S_z|$');
        grid on; xlabel('$t$'); ylabel('$S$'); legend('show');
FMT = '.fig'; save_figure

%% Energy balance
Gn_mid = (G_n(1:end-1)+G_n(2:end))/2.0;
Ga_mid = (G_a(1:end-1)+G_a(2:end))/2.0;
DE_mid = (D_E(1:end-1)+D_E(2:end))/2.0;
DW_mid = (D_W(1:end-1)+D_W(2:end))/2.0;
fig = figure; FIGNAME = ['Energy_balance',sprintf('_%.2d',SID)];
    subplot(221); title('Energy balance terms')
        plot(Ts(2:end),dEdt,'DisplayName','$\partial_t E$'); hold on;
        plot(Ts, G_n, 'DisplayName', '$\Gamma_n$');
        plot(Ts, G_a, 'DisplayName', '$\Gamma_a$');
        plot(Ts, D_E, 'DisplayName', '$D_E$');
        grid on; xlabel('$t$');  legend('show');
    subplot(223); title('Enstrophy balance terms')
        plot(Ts(2:end),dWdt,'DisplayName','$\partial_t W$'); hold on;
        plot(Ts, G_n, 'DisplayName', '$\Gamma_n$');
        plot(Ts, D_W, 'DisplayName', '$D_W$');
        grid on; xlabel('$t$');  legend('show');
    subplot(122); title('Conservation rel. error')
        semilogy(Ts(2:end),100*abs(dEdt - (Gn_mid-Ga_mid-DE_mid))./abs(dEdt),'DisplayName', 'Energy'); hold on;
        plot(Ts(2:end),100*abs(dWdt - (Gn_mid-DW_mid))./abs(dWdt),'DisplayName', 'Enstrophy')
        grid on; xlabel('$t$'); ylabel('$\epsilon[\%]$'); legend('show');
FMT = '.fig'; save_figure

%% Spectra energy
fig = figure; FIGNAME = ['Energy_kin_KZ',sprintf('_%.2d',SID)];
    semilogy(kr(floor(end/2)+1:end),E_kin_KR(floor(end/2)+1:end),'o','DisplayName','$\sum_y\langle|ik\tilde\phi_i|^2\rangle_t$')
    hold on;
    loglog(kz(floor(end/2)+1:end),E_kin_KZ(floor(end/2)+1:end),'o','DisplayName','$\sum_x\langle|ik\tilde\phi_i|^2\rangle_t$')
    grid on; xlabel('$k$');  legend('show');
FMT = '.fig'; save_figure

%% CFL condition
fig = figure; FIGNAME = ['CFL',sprintf('_%.2d',SID)];
    semilogy(Ts,dz./ExB,'-','DisplayName','$|\nabla \phi|\Delta y$');
    hold on;
    plot(Ts,dr*dz/MU*ones(1,numel(Ts)),'-','DisplayName','$\Delta x \Delta y / \mu$');
    plot(Ts,dz/KAPPA*ones(1,numel(Ts)),'-','DisplayName','$\Delta y/ \kappa$');
    plot(Ts,dt*ones(1,numel(Ts)),'--k','DisplayName','$\Delta t$');
    grid on; xlabel('$t$'); ylabel('$\Delta t$'); legend('show');
FMT = '.fig'; save_figure

%% Space-Time diagram at KZ = 0
plt = @(x) real(x);
%% phi
fig = figure; FIGNAME = ['phi_ST',sprintf('_%.2d',SID)];
    [TY,TX] = meshgrid(Ts,z);
    pclr = pcolor(TX,TY,(plt(phi_ST))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$x\,(y=0)$'); ylabel('$t$'); title('$\phi$');
FMT = '.fig'; save_figure

if 0
%% Show frame
it = min(1,numel(Ts));
fig = figure; FIGNAME = ['frame',sprintf('_%.2d',SID)];
    subplot(221); plt = @(x) fftshift((real(x)));
        pclr = pcolor(fftshift(KR),fftshift(KZ),plt(PP(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$\hat\phi$');
    subplot(222); plt = @(x) fftshift(real(x));
        pclr = pcolor(fftshift(KR),fftshift(KZ),plt(NN(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$\hat n$');
    subplot(223); plt = @(x) fftshift(real(x));
        pclr = pcolor(fftshift(KR),fftshift(KZ),plt(ZZ(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$\hat\zeta$');
    subplot(224); plt = @(x) fftshift(log10(abs(x)));
        pclr = pcolor(fftshift(KR),fftshift(KZ),plt(SN(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$');legend('$\hat S_n$');
FMT = '.fig'; save_figure
end
%%
DELAY = 0.067; skip_ = 5;
if 0
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vorticity
GIFNAME = ['zeta',sprintf('_%.2d',SID)]; FIELDNAME = '$\zeta$';
FIELD = real(zeta(:,:,1:skip_:end)); X = XX; Y = YY; T = Ts(1:skip_:end);
create_gif
%% Density
GIFNAME = ['n',sprintf('_%.2d',SID)]; FIELDNAME = '$n$';
FIELD = real(ni(:,:,1:skip_:end)); X = XX; Y = YY; T = Ts(1:skip_:end);
create_gif
%% Phi
GIFNAME = ['phi',sprintf('_%.2d',SID)]; FIELDNAME = '$\phi$';
FIELD = real(phi(:,:,1:skip_:end)); X = XX; Y = YY; T = Ts(1:skip_:end);
create_gif
%% FPhi
GIFNAME = ['Fphi',sprintf('_%.2d',SID)]; FIELDNAME = '$\phi$';
FIELD = real(PP(:,:,1:skip_:end)); X = KR; Y = KZ; T = Ts(1:skip_:end);
create_gif
end