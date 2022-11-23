function t = Imp_opt(Imps, Ls, fm, rho, draw)
    arguments
        Imps (1,:) double
        Ls (1,:) double
        fm (1,1) double = 1
        rho (1,1) double = 1
        draw (1,1) logical = 0
    end

    vin = [1;-1];
    ks = 2*pi*fm*rho./Imps;
    Matrix = [[1,0];[0,1]];
    for i = 1 : size(Ls, 2)
        Matrix = N(ks(i+1),ks(i))*M(ks(i),Ls(i))*Matrix;
    end
    
    r = -Matrix(2,1)/Matrix(2,2);
    t = Matrix(1,1)+Matrix(1,2)*r;
    vin(2) = vin(1)*r;
    t = abs(t);
    if draw
        Imp_draw(Imps, Ls, fm, rho, vin)
    end
end

function matrix = M(k,L)
    %free space propogation
    matrix = [[exp(1i*k*L),0];[0,exp(-1i*k*L)]];
end

function matrix = N(k2,k1)
    %interface
    matrix = 0.5*[[1+k1/k2,1-k1/k2];[1-k1/k2,1+k1/k2]];
end

function Imp_draw(Imps, Ls, fm, rho, vin)
    arguments
        Imps (1,:) double
        Ls (1,:) double
        fm (1,1) double = 1
        rho (1,1) double = 1
        vin (2,1) double = [1,1]
    end
    phase = pi/2;
    Imps(end + 1) = 1;
    ks = 2*pi*fm*rho./Imps;
    n = 100;
    xs = zeros(length(Ls), n+1);
    for i = 1:length(Ls)
        xstart = sum(Ls(1:i)) - Ls(1);
        xend = sum(Ls(1:i));
        xs(i, :) = xstart : (xend - xstart)/n : xend;
    end
    ys = xs;
    Matrix = [[1,0];[0,1]];
    for i = 1:length(Ls)
        vs = Matrix*vin;
        ys(i,:) = exp(1i*ks(i)*(xs(i,:) - xs(i,1)))*vs(1) + exp(-1i*ks(i)*(xs(i,:) - xs(i,1)))*vs(2);
        Matrix = N(ks(i+1),ks(i))*M(ks(i),Ls(i))*Matrix;
    end
    xs = xs';
    ys = ys';
    xs = xs(:);
    ys = real(ys(:)*exp(1i*phase));
    plot(xs, ys);
end