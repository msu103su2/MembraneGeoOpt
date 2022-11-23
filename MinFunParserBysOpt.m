function re = MinFunParserBysOpt(OptObj, x)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    L = 1e-3;
    Z0 = sqrt(OptObj.sigma0*OptObj.rho);
    n = int16(length(table2array(x))/2);
    Imps = table2array(x(1,1:n));
    Imps(end+1) = 1;
    Imps = Imps*Z0;
    Ls = table2array(x(1,n+1:end));
    Ls = Ls*L/sum(Ls);
    OptObj.Imps = Imps;
    OptObj.Ls = Ls;
    
    g = @(fm) (minFun(OptObj, fm));
    [fm, re] = fminbnd(g, 100e3, 10000e3);
    OptObj.Setfm(fm);
end