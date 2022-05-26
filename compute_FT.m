function [Cxx,kr_range] = compute_FT(psi_x,nsigma)

L = size(psi_x,1);

%Compute FT:
A = fft2(psi_x);
A = A/L;
shift = ceil(L/2);
FT = zeros(L,L);
for i = 1 : L
    i_new = mod2(i+shift-1,L);
    for j = 1 : L
        j_new = mod2(j+shift-1,L);
        FT(i_new,j_new) = abs(A(i,j))^2;
    end
end

%Average over wavevectors:
kx_range = (-pi+2*pi()/L):(2*pi()/L):(pi());
ky_range = kx_range;
dkx = (2*pi()/L);
kr_range = dkx:dkx:pi();
sigma = nsigma*dkx;

Cxx = zeros(length(kr_range),1);
for k = 1 : length(kr_range)
    kr = kr_range(k);
    norm = 0;
    for i = 1 : length(kx_range)
        kx = kx_range(i);
        for j = 1 : length(ky_range)
            ky = ky_range(j);
            kmod = sqrt(kx^2+ky^2);
            weight = exp(-(kmod-kr)^2/sigma^2);
            norm = norm+weight;
            Cxx(k) = Cxx(k) + FT(i,j)*weight;
        end
    end
    Cxx(k) = Cxx(k)/norm;
end

end