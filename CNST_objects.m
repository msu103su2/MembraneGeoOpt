classdef CNST_objects
    %CNST_OBJECTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fileID
    end
    
    methods
        function obj = CNST_objects(fileID)
            %CNST_OBJECTS Construct an instance of this class
            %   Detailed explanation goes here
            obj.fileID = fileID;
        end
        
        function Phc_uc(obj, ucGeomParams, xc, yc, resultCell, resultLayer, workingCell,...
            workingLayers, cleanWorkingSpace)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                ucGeomParams (1,1) struct
                xc (1,1) double = 0
                yc (1,1) double = 0
                resultCell (1,1) string = 'top'
                resultLayer (1,:) = 1              
                workingCell (1,1) string = 'top'
                workingLayers (1,:) int8 = [2, 3]
                cleanWorkingSpace (1,1) logical = 0
            end
            
            L = ucGeomParams.L;
            W = ucGeomParams.W;
            r = ucGeomParams.r;
            poleL = ucGeomParams.poleL;
            poleW = ucGeomParams.poleW;
            
            fprintf(obj.fileID, sprintf('<%s struct>\n', workingCell));
            fprintf(obj.fileID, sprintf('%i layer\n', workingLayers(1)));
            fprintf(obj.fileID, sprintf('%.3f %.3f %.3f %.3f 0 rectangleC\n', xc, yc, L, W));
            fprintf(obj.fileID, sprintf('<genArea1 %s %i genArea>\n', workingCell, workingLayers(1)));
            fprintf(obj.fileID, sprintf('%i layer\n', workingLayers(2)));
            fprintf(obj.fileID, sprintf('%.3f %.3f %.3f %.3f %.3f %.3f 0 roundrectC\n', xc, yc, poleL, poleW, r, r));
            fprintf(obj.fileID, sprintf('<genArea2 %s %i genArea>\n', workingCell, workingLayers(2)));
            fprintf(obj.fileID, sprintf('<%s struct>\n', resultCell));
            fprintf(obj.fileID, sprintf('<genArea1 genArea2 %i subtract>\n', resultLayer));
        end
        
        function Phc(obj, phcGeomParams, ucGeomParams_list, xc, yc, resultCell, resultLayer, workingCell,...
            workingLayers, cleanWorkingSpace)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                phcGeomParams (1,1) struct
                ucGeomParams_list (1,:) struct
                xc (1,1) double = 0
                yc (1,1) double = 0
                resultCell (1,1) string = 'top'
                resultLayer (1,:) = 1              
                workingCell (1,1) string = 'top'
                workingLayers (1,:) int8 = [2, 3]
                cleanWorkingSpace (1,1) logical = 0
            end
            dL = phcGeomParams.dL;
            dW = phcGeomParams.dW;
            rTheta = phcGeomParams.rTheta;
            ucXShifts = phcGeomParams.ucXShifts;
            ucYShifts = phcGeomParams.ucYShifts;
            
            fprintf(obj.fileID, sprintf('<%s struct>\n', workingCell));
            fprintf(obj.fileID, sprintf('%i layer\n', workingLayers(1)));
            fprintf(obj.fileID, sprintf('%.3f %.3f %.3f %.3f 0 rectangleC\n', xc, yc, dL, dW));
            
            for i = 1 : length(ucGeomParams_list)
                ucxc = xc + ucXShifts(i);
                ucyc = yc + ucYShifts(i);
                obj.Phc_uc(ucGeomParams_list(i), ucxc, ucyc, workingCell, ...
                    workingLayers(2), strcat(workingCell, '_uc', sprintf('_%iL', i)), workingLayers)
                ucxc = xc - ucXShifts(i);
                ucyc = yc - ucYShifts(i);
                obj.Phc_uc(ucGeomParams_list(i), ucxc, ucyc, workingCell, ...
                    workingLayers(2), strcat(workingCell, '_uc', sprintf('_%iR', i)), workingLayers)
            end
            fprintf(obj.fileID, sprintf('<genArea1 %s %i genArea>\n', workingCell, workingLayers(1)));
            fprintf(obj.fileID, sprintf('<genArea2 %s %i genArea>\n', workingCell, workingLayers(2)));
            fprintf(obj.fileID, sprintf('<%s struct>\n', resultCell));
            fprintf(obj.fileID, sprintf('<genArea1 genArea2 %i or>\n', resultLayer));
        end
    end
end

