F = @(x,y) exp(-x.^2-y.^2);    FF = F(X,Y);
G = @(x,y) sin(x) * exp(y.^2); GG = G(X,Y);

solution = @(c,d) exp(-1/4-d.^2/2)*pi.*sin(c)/sqrt(2);

N = 64;
L = 10;
dx = 2*L/N;
x = linspace(-L,L,N);
y = x;

[Y,X] = meshgrid(y,x);
Pad   = 2.0;
CONV  = conv_thm_2D(FF, GG, Pad)*dx*dx;
%%
CONV_fft2=conv_fft2(FF,GG,'same');
%%
figure; pclr = pcolor(x,y,solution(X,Y)); set(pclr, 'edgecolor','none');
%%
figure; pclr = pcolor(x,y,CONV); set(pclr, 'edgecolor','none');
%%
figure; pclr = pcolor(x,y,CONV_fft2); set(pclr, 'edgecolor','none');
%%
figure; pclr = pcolor(x,y,abs(CONV-solution(X,Y))); set(pclr, 'edgecolor','none');
%%
figure; pclr = pcolor(x,y,(CONV+1e-15)./(solution(X,Y)+1e-15)); set(pclr, 'edgecolor','none');
