function [ A_H ] = hermitian_sym( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[N,N]  = size(A);
M      = 2*N;
A_H    = zeros(M);
A_conj = conj(A);

% I   quadrant
A_H(N+1:M,N+1:M) = A;
% II  quadrant
A_H(  1:N,N+1:M) = flip(A_conj,1);
% III quadrant
A_H(  1:N,  1:N) = flip(flip(A_conj,1),2);
% IV  quadrant
A_H(N+1:M,  1:N) = flip(A_conj,2);

end

