function [params, links] = init_shan()
%import necessary libarary
import com.comsol.model.*
import com.comsol.model.util.*

%----params defined section-----
MS = 3e-6;

%----Model node creation and modification
model = ModelUtil.create('Model');
model.hist.disable;

%----Geom creation----
geom1 = model.geom.create('geom1', 3);

%----geom.wp1---------
wp1 = geom1.feature.create('wp1', 'WorkPlane');

%Prepare a straight string as template

%----Comsol Model params defination
model.param.set('a',[num2str(params.geom.a) '[m]']);
model.param.set('b',[num2str(params.geom.b) '[m]']);
model.param.set('h',[num2str(params.geom.h) '[m]']);
model.param.set('MS',[num2str(params.MS) '[m]']);

%other initializations
exts(1) = geom1.feature.create('ext1','Extrude');
exts(1).selection('input').set('wp1');

Si3N4 = model.material.create('Si3N4');
Si3N4.selection.all;
Si3N4.label('Si3N4');
Si3N4.materialModel('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
Si3N4.materialModel('def').set('thermalexpansioncoefficient', {'2.3e-6[1/K]' '0' '0' '0' '2.3e-6[1/K]' '0' '0' '0' '2.3e-6[1/K]'});
Si3N4.materialModel('def').set('heatcapacity', '700[J/(kg*K)]');
Si3N4.materialModel('def').set('relpermittivity', {'9.7' '0' '0' '0' '9.7' '0' '0' '0' '9.7'});
Si3N4.materialModel('def').set('density', '3100[kg/m^3]');
Si3N4.materialModel('def').set('thermalconductivity', {'20[W/(m*K)]' '0' '0' '0' '20[W/(m*K)]' '0' '0' '0' '20[W/(m*K)]'});
Si3N4.materialModel('def').set('youngsmodulus', '250e9[Pa]');
Si3N4.materialModel('def').set('poissonsratio', '0.23');

mesh = model.mesh.create('mesh', 'geom1');
Msize = mesh.feature('size');
ftri = mesh.feature.create('ftri', 'FreeTri');
swel = mesh.create('swel', 'Sweep');

solid = model.physics.create('solid', 'SolidMechanics', 'geom1');
fix1 = solid.create('fix1', 'Fixed', 2);
iss1 = solid.feature('lemm1').create('iss1', 'InitialStressandStrain', 3);

std = model.study.create('std');
stat = std.feature.create('stat', 'Stationary');
stat.set('geometricNonlinearity', true);
ModelUtil.showProgress(true);
eig = std.feature.create('eig', 'Eigenfrequency');
eig.set('neigsactive', true);
eig.set('geometricNonlinearity', true);
eig.set('useadvanceddisable', true);

links = Links(model, geom1, wp1, exts, mesh, Msize, ftri, swel, swe2, iss1, ...
    fix1, std, eig, solid, ref);
end