sigma = 0.1e9;
rho = 3210;
h = 100e-9;
w = 10e-6;

L0 = 50e-6;
Z0 = sqrt(sigma*rho);

n = 20;
Ls = zeros(1, n) + 0.5;
Imps = zeros(1, n+1) + 1;
Imps(2:2:end) = sqrt(2);
Ls = Imps(1:end-1).*Ls;

fmlo = 100e3;
fmhi = 5000e3;
fms = fmlo:(fmhi - fmlo)/200:fmhi;
ts = zeros(length(fms));
Qs = zeros(length(fms));
ms = zeros(length(fms));
gs = zeros(length(fms));

%Imp_opt(Z0*Imps, L0*Ls, fm, rho)
%find gap frequency;
g = @(fm) (Imp_opt(Z0*Imps, L0*Ls, fm, rho));
[fm, val] = fminbnd(g,100e3, 5000e3);

O1 = ImpOpt(Z0*Imps, L0*Ls, fm);
g = @(fm) (minFun(O1, fm));
[fm, val] = fminbnd(g, 100e3, 5000e3);


for i = 1 : length(fms)
    O1.Setfm(fms(i))
    ts(i) = abs(O1.t);
    ms(i) = O1.Meff();
    Qs(i) = O1.Q_arb();
    gs(i) = g(fms(i));
end
plot(fms, ts);
Imp_opt(Z0*Imps, L0*Ls, fm, rho, 1);

%%%%%%%%%%%%%%%%%
sigma = 1e9;
rho = 3100;
h = 100e-9;
w = 10e-6;

L0 = 100e-6;
Z0 = sqrt(sigma*rho);

n = 20;
Ls = zeros(1, n) + 0.5;
Imps = zeros(1, n+1) + 1;
Imps(2:2:end) = sqrt(2);
Ls = Imps(1:end-1).*Ls;

options = optimoptions('fmincon',...
    'PlotFcn','optimplotfvalconstr',...
    'Display','iter');
options.MaxFunctionEvaluations = 3e4;
O1 = ImpOpt(Z0*Imps, L0*Ls, 1e6);
x0 = [Imps(1:end-1); Ls];

xchache = x0;
minfval = 2e13;
randr = 0.5;
for i = 1:10
    lb = [zeros(1, n) + 0.99; zeros(1, n) + 0.1];
    ub = [zeros(1, n) + sqrt(2.01); zeros(1, n) + 10];
    r = randr*(rand(size(lb))-0.5) + 1;
    %x0  = r .*(ub-lb)+lb;
    xin = x0.*r;
    xin = max(xin, lb*1.01);
    xin = min(xin, ub*0.99);
    g = @(y) (MinFunParser(O1, y));
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    [x,fval] = fmincon(g,xin,A,b,Aeq,beq,lb,ub,nonlcon,options);
    if fval < minfval
        minfval = fval;
        xchache = x;
    end
end

O1.Setvin([1;O1.r]);
O1.ModeShape()
O1.Draw(0)

gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0', xin,...
    'objective',g,'lb',lb,'ub',ub);
x = run(gs,problem);

%%%option1
n =20;
Ls = zeros(1, n) + 0.5;
Imps = zeros(1, n+1) + 1;
Imps(2:2:end) = sqrt(2);
Ls = Imps(1:end-1).*Ls;
NumSeedPoints = 10;
X0s = zeros(NumSeedPoints, 2*n);
lb = [zeros(1, n) + 0.99, zeros(1, n) + 0.1];
ub = [zeros(1, n) + sqrt(2.01), zeros(1, n) + 10];
randr = 0.5;
xin0 = [Imps(1:end-1), Ls];
Y0s = zeros(NumSeedPoints, 1);
for i = 1:NumSeedPoints-1
    r = randr*(rand(size(lb))-0.5) + 1;
    xin = xin0.*r;
    xin = max(xin, lb);
    xin = min(xin, ub);
    Y0s(i) = g(array2table(xin));
    X0s(i,:) = xin;
end
X0s(end,:) = [Imps(1:end-1), Ls];
Y0s(end) = g(array2table(X0s(end,:)));
X0s = array2table(X0s);


%%option2
n =20;
Ls = zeros(1, n) + 0.5;
Imps = zeros(1, n+1) + 1;
Imps(2:2:end) = sqrt(2);
Ls = Imps(1:end-1).*Ls;
NumSeedPoints = 10;
X0s = zeros(NumSeedPoints, 2*n);
lb = [zeros(1, n) + 0.99, zeros(1, n) + 0.1];
ub = [zeros(1, n) + sqrt(2.01), zeros(1, n) + 10];
randr = 0.8;
xin0 = [Imps(1:end-1), Ls];
Y0s = zeros(NumSeedPoints, 1);
for i = 1:NumSeedPoints
    r = randr*(rand(size(lb))-0.5) + 1;
    xin  = r .*(ub-lb)+lb;
    %xin = xin0.*r;
    %xin = max(xin, lb*1.01);
    %xin = min(xin, ub*0.99);
    Y0s(i) = g(array2table(xin));
    X0s(i,:) = xin;
end
%X0s(end,:) = [Imps(1:end-1), Ls];
%Y0s(end) = g(array2table(X0s(end,:)));
X0s = array2table(X0s);

%%option3
n =20;
NumSeedPoints = 100;
X0s = zeros(NumSeedPoints, 2*n);
Y0s = zeros(NumSeedPoints, 1);
for i = 1:NumSeedPoints    
    Ls = zeros(1, n) + 0.5;
    Imps = zeros(1, n+1) + 1;
    Imps(2:2:end) = sqrt(1+rand());
    Ls = Imps(1:end-1).*Ls;
    X0s(i,:) = [Imps(1:end-1), Ls];
    Y0s(i) = g(array2table(X0s(i,:)));
end
X0s = array2table(X0s);

%continue
g = @(y) (MinFunParserBysOpt(O1, y));
BysOptVars = [repmat(optimizableVariable('Imp',[1 sqrt(2)]),1,n), repmat(optimizableVariable('L',[0.1 10]),1,n)];
for i = 1:n
    BysOptVars(1,i) = optimizableVariable(sprintf('Imp%i',i),[1 sqrt(2)]);
    BysOptVars(1,n+i) = optimizableVariable(sprintf('L%i',i),[0.1 10]);
end
x0 = array2table([Imps(1:end-1), Ls]);

results = bayesopt(g, BysOptVars, 'UseParallel', 1, ...
    'NumSeedPoints', NumSeedPoints, 'InitialX', X0s, 'InitialObjective',Y0s, ...
    'MaxObjectiveEvaluations',1200, 'SaveFileName','toy20201020.mat', ...
    'OutputFcn',{@saveToFile});


x = -3:0.01:3;
t = -2:0.01:2;
a= 1;
kx = 2*pi/0.1 + 0*2*pi/a;

uz = cos(2*pi/a*x).^2;
yref = cos(2*pi*t(1) - kx*x).*uz;
pref = plot(x, yref);
hold on
pabs = plot(x, yref.^2);
pref.XDataSource = 'x';
pref.YDataSource = 'yref';
pabs.XDataSource = 'x';
pabs.YDataSource = 'yabs';
for i = 1:length(t)
    yref = cos(2*pi*t(i) - kx*x).*uz;
    yabs = yref.^2;
    refreshdata
    drawnow
end

syms a;
hs = [0:0.01:10];
as = zeros(1, length(hs));
for i = 1:length(hs)
    S = vpasolve(2*a*hs(i)/(hs(i)^2-a^2) == tan(hs(i)), a, [0, inf]);
    as(i) = S;
end

Ls = 182e-9:0.2e-9:185e-9;
spacialOrder = 5;
ris = zeros(spacialOrder*2+1, length(Ls));
tis = zeros(spacialOrder*2+1, length(Ls));
sis = cell(1, length(Ls));
for i = 1 : length(Ls)
    temp = Phc2d(1064e-9, Ls(i), 0/180*pi, 2.6^2, (2.6^2)/2-1, 100e-9, 5);
    [ri,ti,si] = temp.Solve();
    ris(:,i) = ri;
    tis(:,i) = ti;
    sis{i} = si;
end

thetas = -0.5/180*pi:0.1/180*pi:0.5/180*pi;
n = length(thetas);
spacialOrder = 5;
ris = zeros(spacialOrder*2+1, n);
tis = zeros(spacialOrder*2+1, n);
sis = cell(1, n);
for i = 1 : n
    temp = Phc2d(1064e-9, 184.4e-9, thetas(i), 2.6^2, (2.6^2)/2-1, 100e-9, 5);
    [ri,ti,si] = temp.Solve();
    ris(:,i) = ri;
    tis(:,i) = ti;
    sis{i} = si;
end

maxp = max(max(abs(Eg)));
[~,plot] = contourf(temp.meshX, temp.meshZ, real(Eg), 50);
hold on
plot.ZDataSource = 'E';
ts = 0:0.01:1;
for i = 1:length(ts)
    E = real(Eg*exp(1j*2*pi*ts(i)));
    refreshdata
    drawnow
end


test = Phc2d(1064e-9, 653e-9, 1/180*pi, 2.6^2, (2.6^2)/2-1, 100e-9, 10);
test = Phc2d(1064e-9, 640e-9, 1/180*pi, 2.6^2, 0, 100e-9, 10);
test = Phc2d(1064e-9, 653e-9, 1/180*pi, 2.6^2, (2.6^2)/2-1, 100e-9, 10);
test = Phc2d(1064e-9, 653e-9, 0/180*pi, 1, 0, 100e-9, 10);
[ri, ti, si] = test.Solve();
E = test.Exz(ri, ti, si);
%E = flip(E, 1);

test = Phc2d(1064e-9, 653e-9, 2/180*pi, 2.6^2, (2.6^2)/2-1, 100e-9, 10);
test = Phc2d(1064e-9, 653e-9, 0/180*pi, 1, 0, 100e-9, 1);
Eg = test.HGBeamInc(5e-6, 0/180*pi);

cache = Ere;
dataE = cache/max(max(abs(cache)))*256;
fig = figure;
fig.Position = [100 100 700 700];
boxdata = zeros([size(cache),3]);
boxAlphadata = zeros(size(cache));
idx1 = find(test.meshZ == 0);
[~, idx2] = min(abs(test.meshZ-test.d));
boxAlphadata(idx1, :, :) = 255;
boxAlphadata(idx2, :, :) = 255;
img = image(uint8(real(dataE*exp(-1j*2*pi*0))));
hold on;
boximg = image(boxdata);
boximg.AlphaData = boxAlphadata;
colorbar;
tstep = 0.01;
ts = 0:tstep:1;
i = 1;
while true
    i = mod(i, length(ts))+1;
    img.CData = uint8(real(dataE*exp(1j*2*pi*ts(i))));
    drawnow
    pause(tstep)
    if fig.CurrentCharacter > 0
        break;
    end
end

test = Phc2d(1064e-9, 184.4e-9, 1/180*pi, 2.6^2, (2.6^2)/2-1, 100e-9, 10);
[ri, ti, si] = test.Solve();
[Ere, Ein] = test.Exz(ri, ti, si);

Option = 0; %0 for just total without Ein; 1 for Etotal;
cacheIn = Ein;
cache = Ere;
maxE = max([max(max(abs(cache))), max(max(abs(cacheIn)))]);
dataE = cache/maxE*256;
dataEin = cacheIn/maxE*256;
fig = figure;
fig.Position = [100 100 500 500];
boxdata = zeros([size(cache),3]);
boxAlphadata = zeros(size(cache));
idx1 = find(test.meshZ == 0);
[~, idx2] = min(abs(test.meshZ-test.d));
boxAlphadata(idx1, :, :) = 255;
boxAlphadata(idx2, :, :) = 255;
EinAlpha = zeros(size(cache)) + 0.5;
img = image(uint8(real(dataE*exp(1j*2*pi*0))));
hold on;
boximg = image(boxdata);
boximg.AlphaData = boxAlphadata;
imgin = image(uint8(real(dataEin*exp(1j*2*pi*0))));
imgin.AlphaData = EinAlpha;
colorbar;
tstep = 0.01;
ts = 0:tstep:1;
i = 1;
while true
    i = mod(i, length(ts))+1;
    img.CData = uint8(real(dataE*exp(1j*2*pi*ts(i))));
    imgin.CData = uint8(real(dataEin*exp(1j*2*pi*ts(i))));
    drawnow
    pause(tstep)
    if fig.CurrentCharacter > 0
        break;
    end
end

test = Phc2d(1064e-9, 653e-9, 0/180*pi, 2.6^2, 0, 100e-9, 1);
Esic_0 = test.HGBeamInc(5e-6, 0/180*pi);

test = Phc2d(1064e-9, 653e-9, 0.2/180*pi, 2.6^2, 0, 100e-9, 1);
Esic_0_2 = test.HGBeamInc(5e-6, 0.5/180*pi);

test = Phc2d(1064e-9, 653e-9, 0/180*pi, 2.6^2, (2.6^2)/2-1, 100e-9, 1);
Egsic_0 = test.HGBeamInc(5e-6, 0/180*pi);

test = Phc2d(1064e-9, 653e-9, 0.2/180*pi, 2.6^2, (2.6^2)/2-1, 100e-9, 1);
Egsic_0_2 = test.HGBeamInc(5e-6, 0.5/180*pi);


%
x = 0:0.01:5;
tstep = 0.01;
phi = 45/180*pi;
r = 1*exp(-1j*phi);
y = exp(1j*2*pi*x) + r*exp(-1j*2*pi*x);
t = 0:tstep:1;
fig = figure();
p2 = plot(x, y);
ylim([-2,2])
p2.YDataSource = 'yt';
i = 1;
while true
    i = mod(i, length(ts))+1;
    yt = real(y*exp(1j*2*pi*ts(i)));
    refreshdata
    drawnow
    pause(tstep);
    if fig.CurrentCharacter > 0
        break;
    end
end

theta = 1/180*pi;
phc2dobj = test;
x = phc2dobj.meshX;
z = phc2dobj.meshZ;
[X2, Z2] = meshgrid(x, z(z > 0));
if Z2(1,1) < Z2(end,:)
    Z2 = flip(Z2, 1);
    X2 = flip(X2, 1);
end
r = -0.5;
Ek = exp(-1j*k0*(sin(theta)*X2 - cos(theta)*Z2)) + r*exp(-1j*k0*(sin(theta)*X2 + cos(theta)*Z2));

cache = Ek;
%cache = flip(cache, 1);
dataE = cache/max(max(abs(cache)))*256;
fig = figure;
fig.Position = [100 100 500 500];
boxdata = zeros([size(Eg1),3]);
boxAlphadata = zeros(size(Eg1));
idx1 = find(test.meshZ == 0);
[~, idx2] = min(abs(test.meshZ-test.d));
boxAlphadata(idx1, :, :) = 255;
boxAlphadata(idx2, :, :) = 255;
img = image(uint8(real(dataE*exp(-1j*2*pi*0))));
hold on;
boximg = image(boxdata);
boximg.AlphaData = boxAlphadata;
colorbar;
tstep = 0.01;
ts = 0:tstep:1;
i = 1;
while true
    i = mod(i, length(ts))+1;
    img.CData = uint8(real(dataE*exp(1j*2*pi*ts(i))));
    drawnow
    pause(tstep)
    if fig.CurrentCharacter > 0
        break;
    end
end


cache = Ere;
dataE = cache/max(max(abs(cache)))*256;
fig = figure;
fig.Position = [100 100 700 700];
boxdata = zeros([size(cache),3]);
boxAlphadata = zeros(size(cache));
idx1 = find(test.meshZ == 0);
[~, idx2] = min(abs(test.meshZ-test.d));
boxAlphadata(idx1, :, :) = 255;
boxAlphadata(idx2, :, :) = 255;
img = image(uint8(abs(dataE*exp(-1j*2*pi*0))));
hold on;
boximg = image(boxdata);
boximg.AlphaData = boxAlphadata;
colorbar;
tstep = 0.01;
ts = 0:tstep:1;
i = 1;
while true
    i = mod(i, length(ts))+1;
    img.CData = uint8(abs(dataE*exp(1j*2*pi*ts(i))));
    drawnow
    pause(tstep)
    if fig.CurrentCharacter > 0
        break;
    end
end
