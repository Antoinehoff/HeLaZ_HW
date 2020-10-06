function [ conv ] = conv_thm_2D( F, G, Pad )
%conv_thm_2D computes the convolution between F and G with a zero pading Pad
% kspace -> real -> product -> kspace

[Nr, Nz] = size(F);
Mr = Nr * Pad; Mz = Nz * Pad;

f  = ifft2(F,Mr, Mz);
g  = ifft2(G,Mr, Mz);
conv_pad = fft2(f.*g); % convolution becomes product

cmin = (Pad-1)/2; cmax = (Pad+1)/2;

conv = conv_pad(cmin*Nr+1:cmax*Nr, cmin*Nz+1:cmax*Nz); % remove padding
end

