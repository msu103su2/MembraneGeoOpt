function re = MinFunParser(OptObj, x)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    L = 1e-3;
    Z0 = sqrt(OptObj.sigma0*OptObj.rho);
    Imps = x(1,:);
    Imps(end+1) = 1;
    Imps = Imps*Z0;
    Ls = x(2,:);
    Ls = Ls*L/sum(Ls);
    OptObj.Imps = Imps;
    OptObj.Ls = Ls;
    
    g = @(fm) (minFun(OptObj, fm));
    [fm, re] = fminbnd(g, 100e3, 10000e3);
    OptObj.Setfm(fm);
end

