classdef ImpOpt < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Imps
        Ls
        w0 = 10e-6
        h0 = 100e-9
        rho = 3100
        xs
        ys
        ws
        vin = [1;0]
        vout
        sigma0 = 1e9
        Matrix
        r
        t
        fm
        ks
        kb = 1.38e-23;
        T = 300;
        hbar = 1.05e-34;
        rel_step_n = 100;
    end
    
    methods
        function obj = ImpOpt(Imps, Ls, fm)
            obj.Imps = Imps;
            obj.Ls = Ls;
            obj.Setfm(fm);
        end
        
        function rt(obj)         
            obj.r = -obj.Matrix(2,1)/obj.Matrix(2,2);
            obj.t = obj.Matrix(1,1)+obj.Matrix(1,2)*obj.r;
            obj.vout = obj.Matrix*obj.vin;
        end
        
        function ModeShape(obj)
            n = obj.rel_step_n;
            obj.xs = zeros(length(obj.Ls), n+1);
            
            for i = 1:length(obj.Ls)
                if i == 1
                    xstart = 0;
                else
                    xstart = sum(obj.Ls(1:i-1));
                end
                xend = sum(obj.Ls(1:i));
                obj.xs(i, :) = xstart : (xend - xstart)/n : xend;
            end
            obj.ys = obj.xs;
            obj.ws = obj.xs;
            Matrix = [[1,0];[0,1]];
            for i = 1:length(obj.Ls)
                vs = Matrix*obj.vin;
                obj.ys(i,:) = exp(1i*obj.ks(i)*(obj.xs(i,:) - obj.xs(i,1)))*vs(1) + exp(-1i*obj.ks(i)*(obj.xs(i,:) - obj.xs(i,1)))*vs(2);
                Matrix = ImpOpt.N(obj.ks(i+1),obj.ks(i))*ImpOpt.M(obj.ks(i),obj.Ls(i))*Matrix;
                obj.ws(i,:) = obj.w0*obj.Imps(end)^2/obj.Imps(i)^2;
            end
            
            obj.xs = obj.xs';
            obj.ys = obj.ys';
            obj.ws = obj.ws';
            obj.xs = obj.xs(:);
            obj.ys = obj.ys(:);
            obj.ws = obj.ws(:);
        end
        
        function Draw(obj, phase)
            plot(obj.xs, real(obj.ys(:)*exp(1i*phase)), '.');
            hold on
            ws = 1./obj.Imps(1:end-1).^2;
            ws = ws/max(ws);
            for i = 1:length(ws)
                if i == 1
                    xstart = 0;
                else
                    xstart = sum(obj.Ls(1:i-1));
                end
                xend = sum(obj.Ls(1:i));
                xsq = [xstart, xend, xend, xstart, xstart];
                ysq = [-ws(i)/2, -ws(i)/2, ws(i)/2, ws(i)/2, -ws(i)/2];
                plot(xsq, ysq, 'b-', 'LineWidth', 1);
            end
        end

        function Play(obj)
            ws = 1./obj.Imps(1:end-1).^2;
            ws = ws/max(ws);
            fig = figure();
            hold on;
            for i = 1:length(ws)
                if i == 1
                    xstart = 0;
                else
                    xstart = sum(obj.Ls(1:i-1));
                end
                xend = sum(obj.Ls(1:i));
                xsq = [xstart, xend, xend, xstart, xstart];
                ysq = [-ws(i)/2, -ws(i)/2, ws(i)/2, ws(i)/2, -ws(i)/2];
                plot(xsq, ysq, 'b-', 'LineWidth', 1);
            end
            p2 = plot(obj.xs, real(obj.ys(:)*exp(1i*0)), '.');

            ylim([-2,2])
            p2.YDataSource = 'yt';
            i = 1;
            tstep = 0.01;
            ts = 0 : tstep: 1;
            while true
                i = mod(i, length(ts))+1;
                p2.YData = real(obj.ys(:)*exp(1j*2*pi*ts(i)));
                drawnow
                pause(tstep);
                if fig.CurrentCharacter > 0
                    break;
                end
            end
        end
        
        function CalM(obj, fm)
            obj.ks = 2*pi*fm*obj.rho./obj.Imps(1:end-1);
            obj.ks(end+1) = 2*pi*fm*obj.rho*(obj.Imps(end)/obj.Imps(1))^2/obj.Imps(end);
            obj.Matrix = [[1,0];[0,1]];
            for i = 1:length(obj.Ls)
                obj.Matrix = ImpOpt.N(obj.ks(i+1),obj.ks(i))*ImpOpt.M(obj.ks(i),obj.Ls(i))*obj.Matrix;
            end
        end
        
        function Setfm(obj, fm)
            obj.fm = fm;
            obj.CalM(fm);
            obj.ModeShape()
            obj.rt()
        end
        
        function Setvin(obj, vin)
            obj.vin = vin;
            obj.ModeShape()
        end
        
        function m = Meff(obj)
            m = obj.rho*obj.h0*trapz(obj.xs, abs(obj.ys).^2.*obj.ws/(max(abs(obj.ys))^2));
        end
        
        function Q = Q_arb(obj)
            %motion is ys(not normalized), here normalizd by dividing ysmax
            %so U and dU is not directly comparable, but Q is
            ys0 = sin(pi/sum(obj.Ls)*obj.xs);
            m0 = obj.rho*obj.h0*trapz(obj.xs, abs(ys0).^2.*obj.ws/(max(abs(ys0))^2));
            m = obj.Meff();
            U = 2*(pi*obj.fm)^2*m;
            lambda = obj.h0/sum(obj.Ls)*sqrt(80);
            n = obj.rel_step_n;
            eng_clap_ratio = max(abs(obj.ys(end-n+1:end)))/max(abs(obj.ys));
            dU = (mean(obj.ks)^2/sum(obj.Ls))*(eng_clap_ratio^2/lambda*obj.ws(end) + m/m0*sum(obj.Ls)^2*mean(obj.ks)^2/2*mean(obj.ws(end)));
            Q = U/dU;
        end

        function Togds(obj)
            %all in units of um
            cNSTfile = 'string.cnst';
            gdsfile = 'string.gds';
            fileID = fopen(cNSTfile,'w');
            fprintf(fileID,'0.001 gdsReso\n');
            fprintf(fileID,'0.001 shapeReso\n\n');

            xs = [flip(-obj.xs(1:end)'), obj.xs']*1e6;
            ws = [flip(obj.ws(1:end)'), obj.ws']*1e6;
            xsys = zeros(1, 2*length(xs));
            xsys(1:2:end) = xs;
            xsys(2:2:end) = ws/2;
            
            fprintf(fileID,'<top struct>\n');
            fprintf(fileID, sprintf('%i layer\n', 1));
            fprintf(fileID, sprintf('%s 0 0 0 customTaper\n', num2str(xsys)));
            fclose(fileID);
            command = sprintf('java -jar C:\\Users\\shh114\\Documents\\CNSTNanolithographyToolboxV2016.10.01\\CNSTNanolithographyToolboxV2016.10.01.jar cnstscripting %s %s',cNSTfile,gdsfile);
            [status,cmdout] = dos(command)
            interface = Interface(0);
            interface.GDS2DXF(['C:\Users\shh114\','string.gds'], ['Z:\User\Shan\swap\','string.dxf']);
        end
    end
    
    methods(Static)
       function matrix = M(k,L)
            %free space propogation
            matrix = [[exp(1i*k*L),0];[0,exp(-1i*k*L)]];
        end

        function matrix = N(k2,k1)
            %interface
            matrix = 0.5*[[1+k1/k2,1-k1/k2];[1-k1/k2,1+k1/k2]];
        end
    end    
end

