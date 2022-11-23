function meff = Imp_opt_meff(Imps, Ls, w0, h, rho)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    arguments
        Imps (1,:) double
        Ls (1,:) double
        w0 (1,1) double = 10e-6
        h (1,1) double = 100e-3
        rho (1,1) double = 3100
    end
    Ws = w0*sqrt(rho)./Imps(1:end-1).^2;
    meff = 0.5*Ls.*Ws*h*rho;
end

