clear all
%% PARAM AND FUNCTION DEF
F = @(kr,kz) exp((-kr.^2-kz.^2));
G = @(kr,kz) exp((-kr.^2-kz.^2));

L = 10;
N = 64;
kr = linspace(-L,L,N); kz = kr;

[kz, kr] = meshgrid(kr,kz);

FF = F(kr,kz);
GG = G(kr,kz);

xmax = N/(2*(2*L));
xmin = -xmax;
ymax = N/(2*(2*L));
ymin = -ymax;
x = linspace(xmin,xmax,N);
y = linspace(ymin,ymax,N);

%% TEST CONVOLUTION
Pad = 2;
CC = conv_thm_2D(FF,GG, Pad);
figure
subplot(121)
pclr = pcolor(kr,kz,FF);
set(pclr, 'edgecolor','none');
subplot(122)
pclr = pcolor(kr,kz,real(CC));
set(pclr, 'edgecolor','none');

%% TEST IFFT2
Pad = 1;
ff_pad = ifftshift(ifft2(FF,Pad*N, Pad*N));
cmin = (Pad-1)/2; cmax = (Pad+1)/2;
ff = ff_pad(cmin*N+1:cmax*N, cmin*N+1:cmax*N); % remove padding

figure
subplot(121)
pclr = pcolor(kr,kz,FF);
set(pclr, 'edgecolor','none');
subplot(122)
pclr = pcolor(x,y,abs(ff));
set(pclr, 'edgecolor','none');

%% TEST Hermitian symmetry
F = @(kr,kz) kr*exp((-kr.^2-kz.^2));

L = 5;
N = 64;
kr = linspace(0,L,N); kz = kr;

[kz, kr] = meshgrid(kr,kz);

FF = F(kr,kz);
FF_HS = hermitian_sym(FF);

M = 2*N;
Pad = 2;
ff_pad = ifftshift(fft2(FF_HS,Pad*M, Pad*M));
cmin = (Pad-1)/2; cmax = (Pad+1)/2;
ff_sym = ff_pad(cmin*M+1:cmax*M, cmin*M+1:cmax*M); % remove padding
figure
subplot(121)
pclr = pcolor(FF_HS);
set(pclr, 'edgecolor','none');
subplot(122)
pclr = pcolor(imag(ff_sym));
set(pclr, 'edgecolor','none');