classdef Phc2d < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Gamma;
        cutoff_N;
        lambda;
        thetap;
        theta;
        epbar;
        dep;
        d;
        meshX;
        meshZ;
    end

    methods
        function obj = Phc2d(lambda, Gamma, thetap, epbar, dep, d, cutoff_N)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Gamma = Gamma;
            obj.cutoff_N = cutoff_N;
            obj.lambda = lambda;
            obj.thetap = thetap;
            obj.theta = asin(sin(thetap)/sqrt(epbar));
            obj.epbar = epbar;
            obj.dep = dep;
            obj.d = d;
            sized = min([obj.Gamma, obj.d, obj.lambda]);
            span = max([obj.Gamma, obj.d, obj.lambda]);
            span = 5e-6;
            obj.meshX = -80e-6:1e-8:80e-6;
            obj.meshZ = -3e-6:1e-8:3e-6;
        end

        function Eg = HGBeamInc(obj, w0, thetap)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            kx0 = 2*pi/obj.lambda*tan(thetap);
            kcut = 1/(pi*w0);
            ks = -3*kcut:kcut/10:3*kcut;
            ps = zeros(length(ks), 1);
            Eg = zeros(length(obj.meshZ), length(obj.meshX));
            for i = 1:length(ks)
                thetak = atan((ks(i)+ kx0)*obj.lambda/(2*pi));
                ps(i) = exp(-pi^2*w0^2*ks(i)^2);
                obj.setThetap(thetak);
                [ri,ti,si] = obj.Solve();
                Eg = Eg + ps(i)* Exz(obj, ri, ti, si);
                disp(i/length(ks));
            end
        end

        function [re, in] = HGBeamInc2(obj, w0, thetap)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            kx0 = 2*pi/obj.lambda*tan(thetap);
            kcut = 1/(pi*w0);
            kxs = -3*kcut:kcut/10:3*kcut;
            ps = zeros(length(kxs), 1);
            spacialOrder = obj.cutoff_N;
            
            ris = zeros(spacialOrder*2+1, length(ps));
            tis = zeros(spacialOrder*2+1, length(ps));
            sis = cell(1, length(ps));
            cMs = zeros(2*(2*spacialOrder + 1), length(ps));
            thetaps = zeros(length(ps), 1);
            ks = zeros(length(ps), 1);
            lambdaMs = zeros(2*(2*spacialOrder + 1), length(ps));
            w1Ms = zeros(spacialOrder*2+1, 2*(2*spacialOrder + 1), length(ps));
            for i = 1:length(kxs)
                kx = kxs(i)*cos(thetap) + 2*pi/obj.lambda*sin(thetap);
                kz = -kxs(i)*sin(thetap) + 2*pi/obj.lambda*cos(thetap);
                thetak = atan(kx/kz);
                ps(i) = exp(-pi^2*w0^2*kxs(i)^2);
                obj.setThetap(thetak);
                thetaps(i) = thetak;
                [ris(:, i), tis(:, i) ,sis{1, i}, w1Ms(:,:,i), cMs(:, i), lambdaMs(:, i)] = obj.Solve();
                disp(i/length(kxs));
                ks(i) = sqrt(kx^2 + kz^2);
            end
            [re, in] = obj.Exzf(ps, ks, thetaps, ris, tis, cMs, w1Ms, lambdaMs);
        end

        function setThetap(obj, thetap)
            obj.thetap = thetap;
            obj.theta = asin(sin(thetap)/sqrt(obj.epbar));
        end

        function A = MatrixA(obj)
            dimN = 2*obj.cutoff_N + 1;
            m = 2*obj.Gamma*sqrt(obj.epbar)*sin(obj.theta)/obj.lambda;
            orders = -obj.cutoff_N : 1 : obj.cutoff_N;
            a = -2*pi^2*obj.dep/obj.lambda^2;
            bs = 4*pi^2*orders.*(orders - m)/obj.Gamma^2;
            c = 1j*4*pi*(sqrt(obj.epbar)*cos(obj.theta)/obj.lambda);
            M00 = zeros(dimN, dimN);
            M01 = eye(dimN);
            M10 = diag(bs) + diag(repmat(a, 1, dimN-1),1) + diag(repmat(a, 1, dimN-1),-1);
            M11 = eye(dimN) * c;
            A = [M00, M01; M10, M11];
        end

        function [re, in] = Exz(obj, ri, ti, si)
            k1 = 2*pi/obj.lambda;
            k2 = 2*pi/obj.lambda*sqrt(obj.epbar);
            k3 = k1;
            k2z = k2*cos(obj.theta);
            K = 2*pi/obj.Gamma;
            orders = (-obj.cutoff_N :1 : obj.cutoff_N).'; %(n by 1), ri ti si are (1 by n)
            %
            [X2d1, Z2d1] = meshgrid(obj.meshX, obj.meshZ(obj.meshZ < 0));
            [X2d2, Z2d2] = meshgrid(obj.meshX, obj.meshZ(obj.meshZ <= obj.d & obj.meshZ >= 0));
            [X2d3, Z2d3] = meshgrid(obj.meshX, obj.meshZ(obj.meshZ > obj.d));
            %
            if (size(ri,2) == 1)
                ri = ri.';
            end

            if (size(ti,2) == 1)
                ti = ti.';
            end

            k2xi = (k2*sin(obj.theta) - orders*K);
            k1xi = k2xi;
            k3xi = k2xi;
            k1z = sqrt(k1^2 - k1xi.^2);
            k3z = sqrt(k3^2 - k3xi.^2);
            k1z = real(k1z) - 1j*abs(imag(k1z));
            k3z = real(k3z) - 1j*abs(imag(k3z));
            %incmoing
            fin1 = @(x,z) exp(-1j*k1*(sin(obj.thetap)*x + cos(obj.thetap)*z));
            fre1 = @(x,z) ri*exp(-1j*k1xi*x + 1j* k1z*z);
            %f1 = @(x,z) exp(-1j*k1*(sin(obj.thetap)*x + cos(obj.thetap)*z)) + ri*exp(-1j*k1xi*x + 1j*sqrt(k1^2 - k1xi.^2)*z);
            Exzin1 = arrayfun(fin1, X2d1, Z2d1);
            Exzre1 = arrayfun(fre1, X2d1, Z2d1);
            %inside
            fin2 = @(x,z) 0;
            fre2 = @(x,z) (si(z)).'*exp(-1j*k2xi*x - 1j*k2z*z);
            %f2 = @(x,z) (si(z)).'*exp(-1j*k2xi*x - 1j*k2z*z);
            Exzin2 = arrayfun(fin2, X2d2, Z2d2);
            Exzre2 = arrayfun(fre2, X2d2, Z2d2);
            %outgoing
            fin3 = @(x,z) 0;
            fre3 = @(x,z) ti*exp(-1j*k3xi*x - 1j* k3z*(z - obj.d));
            %f3 = @(x,z) ti*exp(-1j*k3xi*x + 1j*sqrt(k3^2 - k3xi.^2)*(z - obj.d));
            Exzin3 = arrayfun(fin3, X2d3, Z2d3);
            Exzre3 = arrayfun(fre3, X2d3, Z2d3);
            re = [Exzre1; Exzre2; Exzre3];
            in = [Exzin1; Exzin2; Exzin3];
        end

        function [r,t,si, w1M, cM, lambdaM] = Solve(obj)
            A = obj.MatrixA();
            delta_order = zeros(2*obj.cutoff_N+1,1);
            delta_order(obj.cutoff_N+1) = 1;
            orders = -obj.cutoff_N : 1 : obj.cutoff_N;
            orders = orders';
            [wM,temp] = eig(A);
            lambdaM = zeros(length(A),1);
            for i = 1:length(A)
                lambdaM(i) = temp(i,i);
            end

            w1M = wM(1:length(wM)/2, :);
            w2M = wM(length(wM)/2+1:end, :);
            k1 = 2*pi/obj.lambda;
            k2 = 2*pi/obj.lambda*sqrt(obj.epbar);
            k3 = k1;
            K = 2*pi/obj.Gamma;

            k2xi = (k2*sin(obj.theta) - orders*K);
            k2z = k2*cos(obj.theta);
            k1xi = k2xi;
            k3xi = k2xi;

            syms RT cM [length(A),1]
            %{
            eqns = [1j*sqrt(k1^2 - k1xi.^2).*(RT(1:length(wM)/2) - delta_order) == w2M*cM - w1M*cM.*(1j*sqrt(k2^2 - k2xi.^2)),...
                RT(length(wM)/2+1:end) == (w1M*(cM.*exp(-lambdaM*obj.d))).*exp(1j*sqrt(k2^2 - k2xi.^2)*obj.d),...
                RT(1:length(wM)/2) + delta_order == w1M*cM,...
                -1j*sqrt(k3^2 - k3xi.^2).*RT(length(wM)/2+1:end) == (w2M*(cM.*exp(-lambdaM*obj.d)) - w1M*(cM.*exp(-lambdaM*obj.d)).*(1j*sqrt(k2^2 - k2xi.^2))).*exp(1j*sqrt(k2^2 - k2xi.^2)*obj.d)];
            
            eqns = [1j*sqrt(k1^2 - k1xi.^2).*(RT(1:length(wM)/2) - delta_order) == w2M*cM +(1j*k2z)*w1M*cM,...
                RT(length(wM)/2+1:end) == (w1M*(cM.*exp(-lambdaM*obj.d)))*exp(-1j*k2z*obj.d),...
                RT(1:length(wM)/2) + delta_order == w1M*cM,...
                -1j*sqrt(k3^2 - k3xi.^2).*RT(length(wM)/2+1:end) == (w2M*(cM.*exp(-lambdaM*obj.d)) + (1j*k2z)*w1M*(cM.*exp(-lambdaM*obj.d)))*exp(-1j*k2z*obj.d)];
            %}
            eqns = [1j*sqrt(k1^2 - k1xi.^2).*(RT(1:length(wM)/2) - delta_order) == w2M*cM-(1j*k2z)*w1M*cM,...
                RT(length(wM)/2+1:end) == (w1M*(cM.*exp(lambdaM*obj.d)))*exp(-1j*k2z*obj.d),...
                RT(1:length(wM)/2) + delta_order == w1M*cM,...
                -1j*sqrt(k3^2 - k3xi.^2).*RT(length(wM)/2+1:end) == (w2M*(cM.*exp(lambdaM*obj.d)) - (1j*k2z)*w1M*(cM.*exp(lambdaM*obj.d)))*exp(-1j*k2z*obj.d)];

            S = solve(eqns, [RT cM]);
            S = structfun(@double,S);
            raysN = obj.cutoff_N*2+1;
            r = S(1:raysN);
            t = S(raysN+1:2*raysN);
            cM = S(2*raysN+1:end);
            si =@(z) w1M*(cM.*exp(lambdaM*z));
        end
        
        function [re, in] = Exzf(obj, ps, ks, thetaps, ris, tis, cMs, w1Ms, lambdaMs)
            thetas = asin(sin(thetaps)/sqrt(obj.epbar));
            k1s = ks;
            k2s = k1s*sqrt(obj.epbar);
            k3s = k1s;
            k2zs = k2s.*cos(thetas);
            K = 2*pi/obj.Gamma;
            
            ps = reshape(ps, [1, length(ps)]);
            thetaps = reshape(thetaps, [1, length(thetaps)]);
            orders = -obj.cutoff_N : 1 : obj.cutoff_N;
            orders = orders';
            
            k2xis = zeros(length(orders), length(ps));
            fsps = cell(1, length(ps));
            for i = 1 : length(ps)
                k2xis(: , i) = k2s(i)*sin(thetas(i)) - orders*K;
            end
            
            k1xis = k2xis;
            k3xis = k2xis;
            
            k1zs = sqrt(repmat(k1s.', [obj.cutoff_N*2+1, 1]).^2 - k1xis.^2);
            k3zs = sqrt(repmat(k3s.', [obj.cutoff_N*2+1, 1]).^2 - k3xis.^2);
            k2zs = repmat(k2zs.', [obj.cutoff_N*2+1, 1]);
            k1zs = real(k1zs) - 1j*abs(imag(k1zs));
            k3zs = real(k3zs) - 1j*abs(imag(k3zs));
            
            fin1 = @(x, z) ps* exp(-1j*(k1s.*( sin(thetaps.')*x + cos(thetaps.')*z )));
            fre1 = @(x, z) ps* (sum(ris.*exp(-1j*k1xis*x + 1j* k1zs*z), 1)).';
            
            fre2 = @(x, z) 0;
            for i = 1 : length(ps)
                fsps{i} =  @(x, z) ps(i) * ((exp(-1j*k2xis(:,i)*x - 1j*k2zs(:,i)*z).')*(w1Ms(:, :, i)*(cMs(:,i).*exp(lambdaMs(:,i)*z))));
                fre2 = @(x, z) fre2(x, z) + fsps{i}(x, z); 
            end
            fin2 = @(x, z) 0;
            
            fin3 = @(x, z) 0;
            fre3 = @(x, z) ps* (sum(tis.*exp(-1j*k3xis*x - 1j* k3zs*(z - obj.d)), 1)).';
            %
            [X2d1, Z2d1] = meshgrid(obj.meshX, obj.meshZ(obj.meshZ < 0));
            [X2d2, Z2d2] = meshgrid(obj.meshX, obj.meshZ(obj.meshZ <= obj.d & obj.meshZ >= 0));
            [X2d3, Z2d3] = meshgrid(obj.meshX, obj.meshZ(obj.meshZ > obj.d));
            %
            
            %incmoing
            Exzin1 = arrayfun(fin1, X2d1, Z2d1);
            Exzre1 = arrayfun(fre1, X2d1, Z2d1);
            %inside
            Exzin2 = arrayfun(fin2, X2d2, Z2d2);
            Exzre2 = arrayfun(fre2, X2d2, Z2d2);
            %outgoing
            Exzin3 = arrayfun(fin3, X2d3, Z2d3);
            Exzre3 = arrayfun(fre3, X2d3, Z2d3);
            re = [Exzre1; Exzre2; Exzre3];
            in = [Exzin1; Exzin2; Exzin3];
        end

        function RayPlot(obj, r, t, si, alpha)
            fig = figure();
            hold on;
            sized = min([obj.Gamma, obj.d, obj.lambda]);
            y = -10*sized:sized/100:10*sized;

            gratingUp = zeros(1, length(y));
            plot([-10*sized, 10*sized], [0,0]);

            gratingUp = zeros(1, length(y))-obj.d;
            plot([-10*sized, 10*sized], [-obj.d,-obj.d]);

            k1 = 2*pi/obj.lambda;
            k2 = 2*pi/obj.lambda*sqrt(obj.epbar);
            k3 = k1;
            K = 2*pi/obj.Gamma;
            
            %plot 0 order
            x = tan(obj.thetap)*y;
            plot(x(y>=0), y(y>=0));

            x = tan(obj.theta)*y;
            plot(x(y<=0 & y>-obj.d), y(y<=0 & y>-obj.d));

            x = tan(obj.thetap)*(y + obj.d);
            plot(x(y<=-obj.d), y(y<=-obj.d));

            if obj.cutoff_N > 0
                for i = 1:obj.cutoff_N
                    p_rip = abs(r(obj.cutoff_N + 1 + i))^2;
                    p_rin = abs(r(obj.cutoff_N + 1 - i))^2;
                    p_tip = abs(t(obj.cutoff_N + 1 + i))^2;
                    p_tin = abs(t(obj.cutoff_N + 1 - i))^2;
                    
                    kx_rip = k2*sin(obj.theta) - i*K;
                    if k1^2 > kx_rip^2
                        kz_rip = sqrt(k1^2 - kx_rip^2);
                        klin_rip = kz_rip/kx_rip;
                        x = y/klin_rip;
                        if alpha
                            plot(x(y>=0), y(y>=0), 'Color',[1, 0, 0, p_rip]);
                        else
                            plot(x(y>=0), y(y>=0), 'Color',[1, 0, 0, 1]);
                        end
                    end
    
                    kx_rin = k2*sin(obj.theta) + i*K;
                    if k1^2 > kx_rin^2
                        kz_rin = sqrt(k1^2 - kx_rin^2);
                        klin_rin = kz_rin/kx_rin;
                        x = y/klin_rin;
                        if alpha
                            plot(x(y>=0), y(y>=0), 'Color',[1, 0, 0, p_rip]);
                        else
                            plot(x(y>=0), y(y>=0), 'Color',[1, 0, 0, 1]);
                        end
                    end
    
                    kx_tip = k2*sin(obj.theta) - i*K;
                    if k3^2 > (kx_rip)^2
                        kz_tip = sqrt(k3^2 - kx_rip^2);
                        klin_tip = kz_tip/kx_tip;
                        x = (y + obj.d)/klin_tip;
                        if alpha
                            plot(x(y<=-obj.d), y(y<=-obj.d), 'Color',[1, 0, 0, p_rip]);
                        else
                            plot(x(y<=-obj.d), y(y<=-obj.d), 'Color',[1, 0, 0, 1]);
                        end
                    end
    
                    kx_tin = k2*sin(obj.theta) + i*K;
                    if k3^2 > (kx_rip)^2
                        kz_tin = sqrt(k3^2 - kx_tin^2);
                        klin_tin = kz_tin/kx_tin;
                        x = (y + obj.d)/klin_tin;
                        if alpha
                            plot(x(y<=-obj.d), y(y<=-obj.d), 'Color',[1, 0, 0, p_rip]);
                        else
                            plot(x(y<=-obj.d), y(y<=-obj.d), 'Color',[1, 0, 0, 1]);
                        end
                    end
                end
            end
            hold off;
        end
    end
end