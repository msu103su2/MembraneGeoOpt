function r = test(k1, k2, phi)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    r = (1-k1)/(1+k1) - (2+4*k1)/(1-k1^2)*(1/(1-(1-k1)/(1+k1)*(k1-k2)/(k1+k2)*exp(1i*phi))-1);
end