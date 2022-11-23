function [gdsfile] = CNST_arm(params, interface)
    ucGeomParams_list = params.ucGeomParams_list;
    phcGeomParams = params.phcGeomParams;
    %all in units of um
    cNSTfile = 'arm.cnst';
    gdsfile = 'arm.gds';
    fileID = fopen(cNSTfile,'w');
    fprintf(fileID,'0.001 gdsReso\n');
    fprintf(fileID,'0.001 shapeReso\n\n');
    
    fprintf(fileID,'<top struct>\n');
    CNST_writter = CNST_objects(fileID);
    CNST_writter.Phc(phcGeomParams, ucGeomParams_list, 0, 0, 'top', 1, 'phc');
    command = sprintf('java -jar C:\\Users\\shh114\\Documents\\CNSTNanolithographyToolboxV2016.10.01\\CNSTNanolithographyToolboxV2016.10.01.jar cnstscripting %s %s',cNSTfile,gdsfile);
    [status,cmdout] = dos(command)
    interface.GDS2DXF(['C:\Users\shh114\','arm.gds'], ['Z:\User\Shan\swap\','arm.dxf']);
end