%% Compute the Si00 term with matlab methods and compare it with HeLaZ
% This script is launched apr�s a special HeLaZ linear run where Si00 is
% still computed.
%% Compute vorticity non linear term {phi,zeta}
%padding
Pad_M = 2.0;
Sz_M = zeros(size(ni));
F_   = zeros(numel(kr), numel(kz));
G_   = zeros(numel(kr), numel(kz));
[KZ, KR] = meshgrid(kz,kr);
for it = 1:numel(time)
    % first term
    F_ = 1i*KR.*phi(:,:,it);
    %Second conv term
    G_ = 1i*KZ.*F_.^2;%1i*KZ.*ze(:,:,it);

    %Conv theorem
    Sz_M(:,:,it) = conv_thm_2D(F_,G_,Pad_M);
    
    % second term
    %First conv term
    F_ = 1i*KZ.*phi(:,:,it);
    %Second conv term
    G_ = 1i*KR.*F_.^2;%1i*KR.*ze(:,:,it);
    %Conv theorem
    Sz_M(:,:,it) = Sz_M(:,:,it) - conv_thm_2D(F_,G_,Pad_M);
end
%% Compute Density non linear term {phi,n}
%padding
Pad_M = 2.0;
Sn_M = zeros(size(ni));
F_   = zeros(numel(kr), numel(kz));
G_   = zeros(numel(kr), numel(kz));
[KZ, KR] = meshgrid(kz,kr);
for it = 1:numel(time)
    % first term
    F_ = 1i*KR.*phi(:,:,it);
    %Second conv term
    G_ = 1i*KZ.*ni(:,:,it);

    %Conv theorem
    Sn_M(:,:,it) = conv_thm_2D(F_,G_,Pad_M);
    
    % second term
    %First conv term
    F_ = 1i*KZ.*phi(:,:,it);
    %Second conv term
    G_ = 1i*KR.*ni(:,:,it);
    %Conv theorem
    Sn_M(:,:,it) = Sn_M(:,:,it) - conv_thm_2D(F_,G_,Pad_M);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison with HeLaZ
Sipj        = load_5D_data(filename, 'Sipj');
Sn_H    = squeeze(Sipj(1,1,:,:,:));
M_mean  = mean(mean(mean(Sn_M)));
H_mean  = mean(mean(mean(Sn_H)));
disp(['Mean^3 Matlab - HeLaz Sn : ',num2str(M_mean-H_mean)])

Sepj        = load_5D_data(filename, 'Sipj');
Sz_H    = squeeze(Sipj(1,1,:,:,:));
M_mean  = mean(mean(mean(Sz_M)));
H_mean  = mean(mean(mean(Sz_H)));
disp(['Mean^3 Matlab - HeLaz Sz : ',num2str(M_mean-H_mean)])
%% Show first and last frame
S_H  = Sz_H;
S_M  = Sz_M;
spec = 'z';
it   = 2;
fig = figure;
    subplot(221)
        pclr = pcolor(kr,kz,log10(abs(S_H(:,:,it))));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('kr');ylabel('kz');
        title(['$log|S_{HeLaZ}^{',spec,'}|$, $t \approx$', sprintf('%.3d',ceil(time(it)))])
    subplot(222)
        pclr = pcolor(kr,kz,log10(abs(S_M(:,:,it))));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('kr');ylabel('kz');
        title(['$log|S_{Matlab}^{',spec,'}|$, $t \approx$', sprintf('%.3d',ceil(time(it)))])
    subplot(223)
        pclr = pcolor(kr,kz,log10(abs(S_H(:,:,end))));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('kr');ylabel('kz');
        title(['$log|S_{HeLaZ}^{',spec,'}|$, $t \approx$', sprintf('%.3d',ceil(time(end)))])
    subplot(224)
        pclr = pcolor(kr,kz,log10(abs(S_M(:,:,end))));
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('x');ylabel('y');
        title(['$log|S_{Matlab}^{',spec,'}|$, $t \approx$', sprintf('%.3d',ceil(time(end)))])
FIGNAME = 'check_Si00';

%% Error first and last frame
ERR_str = log10(abs(  S_H(:,:,1) -   S_M(:,:,1)));
ERR_end = log10(abs(S_H(:,:,end) - S_M(:,:,end)));
fig = figure;
     subplot(211)
        pclr = pcolor(kr,kz,ERR_str);
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('kr');ylabel('kz');
        title(['$log|S_{HeLaZ}^{',spec,'}-S_{Matlab}^{',spec,'}|$, $t \approx$', sprintf('%.3d',ceil(time(1)))])
     subplot(212)
        pclr = pcolor(kr,kz,ERR_end);
        set(pclr, 'edgecolor','none');
        colormap jet; colorbar;
        xlabel('kr');ylabel('kz');
        title(['$log|S_{HeLaZ}^{',spec,'}-S_{Matlab}^{',spec,'}|$, $t \approx$', sprintf('%.3d',ceil(time(end)))])