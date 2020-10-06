function [ A_filtered ] = epsilon_filter( A, epsilon )
%epsilon_filter Filter every small value (<epsilon) and replace by 0

A_filtered = A;
for i = 1,numel(A);
    if abs(A(i)) < epsilon
        A_filtered(i) = 0.0;
    else
        A_filtered(i) = A(i);
    end
    
end

