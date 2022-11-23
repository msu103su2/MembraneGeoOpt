temp = Interface(1);
dxfDir = 'Z:\User\Shan\swap\';
combined = CNST_tilt();
pad = CNST_pad();
arm = CNST_arm();

temp.GDS2DXF(['C:\Users\shh114\',combined], ['Z:\User\Shan\swap\','combined.dxf'])
temp.GDS2DXF(['C:\Users\shh114\',pad], ['Z:\User\Shan\swap\','pad.dxf'])
temp.GDS2DXF(['C:\Users\shh114\',arm], ['Z:\User\Shan\swap\','arm.dxf'])

GeomParams.L = 15;
GeomParams.W = 2;
GeomParams.r = 0.2;
GeomParams.poleL = 12;
GeomParams.poleW = 1;
numofuc = 15;
ucGeomParams_list = repmat(GeomParams, 1, numofuc);
phcGeomParams.dL = 25;
phcGeomParams.dW = 2;
phcGeomParams.ucXShifts = (phcGeomParams.dL+GeomParams.L)/2:GeomParams.L:(phcGeomParams.dL+GeomParams.L)/2 + (numofuc-1)*30;
phcGeomParams.ucYShifts = zeros(1, numofuc) + 0;
phcGeomParams.rTheta = 0;

params.ucGeomParams_list = ucGeomParams_list;
params.phcGeomParams = phcGeomParams;