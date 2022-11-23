classdef Interface
    %INTERFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Model
    end
    
    methods
        function obj = Interface(Model)
            %INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Model = Model;
        end
        
        function GDS2DXF(obj, GDSfile, DXFfile)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            cmd = sprintf("python GDS2DXF.py %s %s", GDSfile, DXFfile);
            system(cmd);
        end
    end
end

