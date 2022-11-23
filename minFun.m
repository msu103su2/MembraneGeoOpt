function cost = minFun(OptObj, fm)
    OptObj.Setfm(fm);
    OptObj.Setvin([1,;OptObj.r]);
    cost = abs(OptObj.t);
    %cost = abs(OptObj.t)*2*pi*fm/OptObj.Q_arb;
end