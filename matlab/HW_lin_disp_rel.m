%% Hasegawa-Wakatani linear dispersion relation
function [maxG, GG] = HW_lin_disp_rel(alpha,kappa,mu,KX,KY)
    gamma = @(kx,ky) max(imag(roots([...
        (kx^2+ky^2),...
        1i*(alpha*((kx^2+ky^2)+1) + mu*(kx^2+ky^2)^2- mu*(kx^2+ky^2)),...
        mu^2*(kx^2+ky^2)^2 - 1i*alpha*ky*kappa])));
%     gamma = @(kx,ky) max(imag(roots([(kx^2+ky^2), 1i*alpha*((kx^2+ky^2)+1), -alpha^2 + alpha - 1i*ky*kappa])));

    %% Frequency space %%%%%%%%%%%%%%%%%%%%%%%%%%
    [Nx_,Ny_] = size(KX);

    %% Evaluation of gamma on the grid %%%%%%%%%%
    GG = zeros(Nx_,Ny_);

    if alpha == 0 || kappa == 0
        maxG = 0;
    else
        for ikx = 1:Nx_
            for iky = 1:Ny_
                kx_ = KX(ikx,iky); ky_ = KY(ikx,iky);
                GG(ikx,iky) = gamma(kx_,ky_);
            end
        end
        maxG = (max(max(GG)));
    end
end