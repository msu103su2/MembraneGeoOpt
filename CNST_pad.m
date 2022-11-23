function [gdsfile] = CNST_pad()
    %all in units of um
    cNSTfile = 'pad.cnst';
    gdsfile = 'pad.gds';
    fileID = fopen(cNSTfile,'w');
    fprintf(fileID,'0.001 gdsReso\n');
    fprintf(fileID,'0.001 shapeReso\n\n');
    
    tiltCenterL = 30;
    
    clapL = 5;
    
    X = -tiltCenterL/2:tiltCenterL/100:tiltCenterL/2;
    Y = 0*(1+cos(3*pi*X/tiltCenterL)) + 5;
    
    %{
    clapXL = -tiltCenterL/2 - clapL : clapL/50 : -tiltCenterL/2;
    clapXR = tiltCenterL/2 : clapL/50 : tiltCenterL/2 + clapL;
    clapYL = clapXL - clapXL(1);
    clapYR =  - (clapXR - clapXR(end));
    
    
    
    X = [clapXL, X, clapXR];
    Y = [clapYL, Y, clapYR];
    %}
    temp = zeros(1, 2*length(X));
    temp(1:2:end) = X;
    temp(2:2:end) = Y;
    coordXY = num2str(temp);
    
    fprintf(fileID,'# Creating tilt part\n');
    fprintf(fileID,'<top struct>\n');
    fprintf(fileID,'1 layer\n');
    fprintf(fileID,sprintf('%s 0 0 0 customTaper\n', coordXY));
    
    command = sprintf('java -jar C:\\Users\\shh114\\Documents\\CNSTNanolithographyToolboxV2016.10.01\\CNSTNanolithographyToolboxV2016.10.01.jar cnstscripting %s %s',cNSTfile,gdsfile);
    [status,cmdout] = dos(command)
end