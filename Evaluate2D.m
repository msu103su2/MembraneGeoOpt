classdef Evaluate2D < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        w;%displacement along z, outofplane motion
        PoissonR;
        PoissonM;
        E;
        MeshX;
        MeshY;
    end
    
    methods
        function obj = Evaluate2D()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function Q = Cal_Q(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Q = 1;
        end
        
        function Set_Possion(obj, PossionR)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.PossionR = PossionR;
            obj.PoissonM = [1,          PossionR,       0;              ...
                            PossionR,   1,              0;              ...
                            0,          0,              (1-PossionR)/2] ...
                            /(1 - PossionR^2);
        end
        
        function Set_YoungModu(obj, YoungModu)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.E = YoungModu;
        end
        
        function Qs = Cal_Q_comsol(obj, model)
            %METHOD1 Summary of this method goes here
            %   Current version has avg Q about 4/3 times the comsol calculation for thin rectangular plate 
            expr_dU = 'z^2*pi*imag(solid.E)/(1-solid.nu^2)*((wXX + wYY)^2 - 2*(1 - solid.nu)*(wXX*wYY - wXY^2))';%center at z=0
            %expr_dU = 'pi*imag(solid.E)*h^3/(12*(1-solid.nu^2))*(wXX + wYY)^2';
            expr_U = '2*solid.rho*pi^2*freq^2*w^2';
            exprs = {expr_dU, expr_U};
            [dUs, Us] = mphint2(model, exprs, 3);
            Qs = 2*pi*Us./dUs;
        end
    end
end

