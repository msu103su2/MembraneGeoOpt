function SingleRunResult = reconstruct(params, links)
import com.comsol.model.*
import com.comsol.model.util.*
eps = 1e-10;

% reset defect params and others
links.model.param.set('a',[num2str(params.geom.a) '[m]']);
links.model.param.set('b',[num2str(params.geom.b) '[m]']);
links.model.param.set('h',[num2str(params.geom.h) '[m]']);
links.model.param.set('MS',[num2str(params.MS) '[m]']);
links.model.material('Si3N4').materialModel('def').set('youngsmodulus', [num2str(params.youngsModulus),'[Pa]']);

%recreate workplane, recreate defect
for i = 1:size(links.exts)
    links.geom1.feature.remove(links.exts(i).tag);
end
links.exts = [];
links.geom1.feature.remove('wp1');
links.wp1 = links.geom1.feature.create('wp1', 'WorkPlane');
links.geom1.feature.move('wp1',0);
links.wp1.set('quickplane', 'xy');

%import dxf
imp1 = links.wp1.geom.create('imp1', 'Import');
imp1.set('type', 'dxf');
imp1.set('filename', DXFfile);
imp1.set('alllayers', {'L1D0' 'L11D0' 'L0D0' 'L2D0'});
imp1.set('layerselection', 'selected');
imp1.set('layers', {'L11D0'});
imp1.set('repairtoltype', 'relative');
imp1.set('repairtol', 1.0E-9);
sca1 = links.wp1.geom.create('sca1', 'Scale');
sca1.set('factor', '1e-6');
sca1.selection('input').set({'imp1'});
%extrude
%it seems from comsol, for 2D object on wp1, there only being 1 face with
%indice 1

links.geom1.run;

%prestress condiciton
links.iss1.set('Sil', cellstr(string(reshape(params.stressTensor,1,[]))));

%study node
links.eig.set('neigs', 100);

%boundary condition
bnd1box = [leftEndCoord(1)-eps, leftEndCoord(1)+eps;...
    -params.defect.B.width*5, params.defect.B.width*5;...
    -eps, params.defect.B.height*5];
bnd2box = [rightEndCoord(1)-eps, rightEndCoord(1)+eps;...
    -params.defect.B.width*5, params.defect.B.width*5;...
    -eps, params.defect.B.height*5];
ftribox = [leftEndCoord(1)-eps, rightEndCoord(1)+eps;...
    -params.defect.B.width*5, params.defect.B.width*5;...
    -eps, eps];
%{
sweldestibox = [leftEndCoord(1)-eps, rightEndCoord(1)+eps;...
    -params.defect.B.width*5, params.defect.B.width*5;...
    -eps, eps];
    %}
idx_bnd1 = mphselectbox(links.model,'geom1', bnd1box, 'boundary');
idx_bnd2 = mphselectbox(links.model,'geom1', bnd2box, 'boundary');
idx_ftri = mphselectbox(links.model,'geom1', ftribox, 'boundary');
%idx_sweldestiface = mphselectbox(links.model,'geom1', sweldestibox, 'boundary');
links.fix1.selection.set([idx_bnd1 idx_bnd2]);
links.Msize.set('hmax', 'MS');
links.Msize.set('hmin', 'MS/4');
links.Msize.set('hcurve', '0.2');
links.ftri.selection.set(idx_ftri);

idx_ext1 = mphselectbox(links.model,'geom1', ftribox+[0,0;0,0;0,heights(1)], 'domain');
links.swel.selection.geom('geom1', 3);
links.swel.selection.set(idx_ext1);
links.swel.create('dis1', 'Distribution');
links.swel.feature('dis1').set('numelem', 2);

if (size(heights, 2) == 2)
    idx_ext2 = mphselectbox(links.model,'geom1', ftribox-[0,0;0,0;heights(2),0], 'domain');
    links.swe2.selection.geom('geom1', 3);
    links.swe2.selection.set(idx_ext2);
    links.swe2.create('dis1', 'Distribution');
    links.swe2.feature('dis1').set('numelem', 2);
end

links.mesh.run;
links.std.run;
[localmodefreq, localmodeEffMass] = Localmode_center_tilt(links,params,lowercenterline);

%export data
SingleRunResult = evaluategeom(links,params,localmodefreq,localmodeEffMass,leftEndCoord, rightEndCoord);
SingleRunResult.EffMass = localmodeEffMass;
SingleRunResult.modes = sortModes(links, params, lowercenterline);
links.swel.feature.remove('dis1');
end