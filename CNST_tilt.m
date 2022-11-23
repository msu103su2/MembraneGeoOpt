function [gdsfile] = CNST_tilt(params, interface)
    %all in units of um
    cNSTfile = 'combined.cnst';
    gdsfile = 'combined.gds';
    fileID = fopen(cNSTfile,'w');
    fprintf(fileID,'0.001 gdsReso\n');
    fprintf(fileID,'0.001 shapeReso\n\n');
    
    ucGeomParams_list = params.ucGeomParams_list;
    phcGeomParams = params.phcGeomParams;
    %all in units of um
    xc = 0;
    yc = 0;
    
    armAngle = 60;
    
    tiltCenterL = 10;
    tiltCenterW = 10;
    X = -tiltCenterL/2:tiltCenterL/100:tiltCenterL/2;
    Y = 0*(1+cos(3*pi*X/tiltCenterL)) + tiltCenterW/2;
    temp = zeros(1, 2*length(X));
    temp(1:2:end) = X;
    temp(2:2:end) = Y;
    coordXY = num2str(temp);
    
    fprintf(fileID,'<top struct>\n');
    fprintf(fileID,'2 layer\n');
    fprintf(fileID,sprintf('%s 0 0 0 customTaper\n', coordXY));
    CNST_writter = CNST_objects(fileID);
    CNST_writter.Phc(phcGeomParams, ucGeomParams_list, 0, 0, 'phc', 1, 'phc_ws');
    
    fprintf(fileID,'<phc struct>\n');
    fprintf(fileID,'2 layer\n');
    fprintf(fileID, sprintf('%.3f %.3f %.3f %.3f 0 rectangleLH\n', xc, yc - phcGeomParams.dW*5, phcGeomParams.ucXShifts(end)*2, phcGeomParams.dW*10));
    fprintf(fileID, sprintf('<genArea1 phc 1 genArea>\n'));
    fprintf(fileID, sprintf('<genArea2 phc 2 genArea>\n'));
    fprintf(fileID, sprintf('<%s struct>\n', 'phc_ins'));
    fprintf(fileID, sprintf('<genArea1 genArea2 %i and>\n', 1));
    
    fprintf(fileID, sprintf('<%s struct>\n', 'top'));
    fprintf(fileID,'3 layer\n');
    fprintf(fileID, sprintf('<phc_ins %.3f %.3f N 1 %.3f instance>\n', X(1) + phcGeomParams.dW*sin(armAngle/180*pi)/2, Y(1) - phcGeomParams.dW*cos(armAngle/180*pi)/2, 180-armAngle));
    fprintf(fileID, sprintf('<phc_ins %.3f %.3f N 1 %.3f instance>\n', X(1) + phcGeomParams.dW*sin(armAngle/180*pi)/2, -Y(1) + phcGeomParams.dW*cos(armAngle/180*pi)/2, 180+armAngle));
    fprintf(fileID, sprintf('<phc_ins %.3f %.3f N 1 %.3f instance>\n', X(end) - phcGeomParams.dW*sin(armAngle/180*pi)/2, Y(end) - phcGeomParams.dW*cos(armAngle/180*pi)/2, armAngle));
    fprintf(fileID, sprintf('<phc_ins %.3f %.3f N 1 %.3f instance>\n', X(end) - phcGeomParams.dW*sin(armAngle/180*pi)/2, -Y(end) + phcGeomParams.dW*cos(armAngle/180*pi)/2, 360-armAngle));
    fprintf(fileID, sprintf('<genArea1 top 2 genArea>\n'));
    fprintf(fileID, sprintf('<genArea2 top 3 genArea>\n'));
    fprintf(fileID, sprintf('<genArea1 genArea2 %i or>\n', 1));
    
    command = sprintf('java -jar C:\\Users\\shh114\\Documents\\CNSTNanolithographyToolboxV2016.10.01\\CNSTNanolithographyToolboxV2016.10.01.jar cnstscripting %s %s',cNSTfile,gdsfile);
    [status,cmdout] = dos(command)
    interface.GDS2DXF(['C:\Users\shh114\', gdsfile], ['Z:\User\Shan\swap\','device.dxf']);
end