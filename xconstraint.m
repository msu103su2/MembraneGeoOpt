function tf = xconstraint(X, ub, lb)
    tf = X<=ub & X>=lb;
end