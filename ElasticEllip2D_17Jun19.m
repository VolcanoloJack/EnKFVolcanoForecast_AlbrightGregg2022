function out = model
%
% ElasticEllip2D_17Jun19.m
%
% Model exported on Oct 11 2022, 14:07 by COMSOL 6.0.0.318.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/patriciagregg/Box Sync/MatlabCode/DataPaper2022/Assimilation');

model.label('2DElasticEllip_17Jun19.mph');

model.param.set('r_dist', '50 [km]', 'model space radius');
model.param.set('z_dist', '50 [km]', 'model space depth');
model.param.set('r1', '1000 [m]', 'reservoir half height');
model.param.set('r2', '2000 [m]', 'reservoir half width');
model.param.set('dP', '100 [MPa]', 'reservoir overpressure');
model.param.set('dZ', '4 [km]', 'depth to reservoir center');
model.param.set('E', '75 [GPa]', 'Youngs modulus');
model.param.set('nu', '0.25', 'Poissons ratio');
model.param.set('G', 'E/(2*(1+nu))', 'shear modulus');
model.param.set('K', 'E/(3*(1-2*nu))', 'bulk modulus');
model.param.set('rhor', '2500 [kg/m^3]', 'host rock density');
model.param.set('g', '9.81 [m/s^2]', 'gravity');
model.param.set('aif', '25 [deg]', 'MC angle of internal friction');

model.component.create('comp1', false);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').curvedInterior(false);

model.component('comp1').geom('geom1').axisymmetric(true);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').repairTolType('relative');
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', {'0' '-z_dist'});
model.component('comp1').geom('geom1').feature('r1').set('size', {'r_dist' 'z_dist'});
model.component('comp1').geom('geom1').create('e1', 'Ellipse');
model.component('comp1').geom('geom1').feature('e1').set('pos', {'0' '-dZ'});
model.component('comp1').geom('geom1').feature('e1').set('semiaxes', {'r2' 'r1'});
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'0' '-dZ-r1-500'});
model.component('comp1').geom('geom1').feature('r2').set('size', {'2*r2' 'dZ+r1+500'});
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r1' 'r2'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'e1'});
model.component('comp1').geom('geom1').create('ls1', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls1').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls1').set('coord2', [25000 0]);
model.component('comp1').geom('geom1').feature('ls1').selection('vertex1').set('dif1(1)', 5);
model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'relative');
model.component('comp1').geom('geom1').run;

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('MagmaLd', 'dP + rhor*g*-z', 'Magma load on reservoir wall');
model.component('comp1').variable('var1').set('beta', '(z + dZ) / r1');
model.component('comp1').variable('var1').set('alpha', 'atan((z+dZ)/r)');
model.component('comp1').variable('var1').set('tau', '(solid.sp3-solid.sp1)/2 * cos(aif)', 'shear stress of MC plane');
model.component('comp1').variable('var1').set('sigma', '(solid.sp3+solid.sp1)/2 - (solid.sp3-solid.sp1)/2 * sin(aif)', 'normal stress of MC plane');
model.component('comp1').variable('var1').set('MCFail', 'abs(tau) - (abs(sigma) * tan(aif))', 'MC failure if > Cohesion');
model.component('comp1').variable('var1').set('sigmaX', '(cos((aif+pi/2)/2))^2*solid.sp3+(cos((aif-pi/2)/2)^2)*solid.sp1', 'long for sigma');
model.component('comp1').variable('var1').set('tauX', 'sqrt(((cos((aif-pi/2)/2)*cos((aif+pi/2)/2))^2)*(solid.sp3-solid.sp1)^2)', 'long form tau');
model.component('comp1').variable('var1').set('MCFailX', 'abs(tauX) - (abs(sigmaX) * tan(aif))');

model.view.create('view2', 3);
model.view.create('view3', 3);

model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');
model.component('comp1').physics('solid').feature('lemm1').create('iss1', 'InitialStressandStrain', 2);
model.component('comp1').physics('solid').create('roll1', 'Roller', 1);
model.component('comp1').physics('solid').feature('roll1').selection.set([2 10]);
model.component('comp1').physics('solid').create('bndl1', 'BoundaryLoad', 1);
model.component('comp1').physics('solid').feature('bndl1').selection.set([11 12]);
model.component('comp1').physics('solid').create('gr1', 'Gravity', 2);
model.component('comp1').physics('solid').feature('gr1').selection.set([1 2]);

model.component('comp1').mesh('mesh1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').create('size2', 'Size');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 1);
model.component('comp1').mesh('mesh1').feature('size1').selection.set([11 12]);
model.component('comp1').mesh('mesh1').feature('size2').selection.geom('geom1', 1);
model.component('comp1').mesh('mesh1').feature('size2').selection.set([5 6 8]);

model.thermodynamics.label('Thermodynamics Package');

model.component('comp1').view('view1').axis.set('xmin', -54044.25);
model.component('comp1').view('view1').axis.set('xmax', 104044.2421875);
model.component('comp1').view('view1').axis.set('ymin', -68951.328125);
model.component('comp1').view('view1').axis.set('ymax', 18951.32421875);

model.component('comp1').physics('solid').prop('ShapeProperty').set('order_displacement', 2);
model.component('comp1').physics('solid').feature('lemm1').set('IsotropicOption', 'KG');
model.component('comp1').physics('solid').feature('lemm1').set('E_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('K_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('K', 'K');
model.component('comp1').physics('solid').feature('lemm1').set('G_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('G', 'G');
model.component('comp1').physics('solid').feature('lemm1').set('rho_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('rho', 'rhor');
model.component('comp1').physics('solid').feature('lemm1').set('minput_temperature_src', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('minput_strainreferencetemperature_src', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').feature('iss1').set('Sil', {'rhor*g*z'; '0'; '0'; '0'; 'rhor*g*z'; '0'; '0'; '0'; 'rhor*g*z'});
model.component('comp1').physics('solid').feature('dcnt1').set('pairDisconnect', true);
model.component('comp1').physics('solid').feature('dcnt1').label('Contact');
model.component('comp1').physics('solid').feature('dcont1').set('pairDisconnect', true);
model.component('comp1').physics('solid').feature('dcont1').label('Continuity');
model.component('comp1').physics('solid').feature('roll1').set('NormalAxis', 'zAxis');
model.component('comp1').physics('solid').feature('bndl1').set('LoadType', 'FollowerPressure');
model.component('comp1').physics('solid').feature('bndl1').set('FollowerPressure', 'MagmaLd');
model.component('comp1').physics('solid').feature('gr1').set('g', {'0'; '0'; '-g'});

model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size1').set('hmax', 75);
model.component('comp1').mesh('mesh1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('size1').set('hgrad', 1.5);
model.component('comp1').mesh('mesh1').feature('size1').set('hgradactive', false);
model.component('comp1').mesh('mesh1').feature('size2').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size2').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size2').set('hmax', 300);
model.component('comp1').mesh('mesh1').feature('size2').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').set('smoothmaxiter', 8);
model.component('comp1').mesh('mesh1').feature('ftri1').set('smoothmaxdepth', 8);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('param', 'Parametric');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').label('Parametric Solutions 1');
model.sol.create('sol6');

model.result.dataset.remove('dset3');
model.result.dataset('dset2').set('solution', 'sol6');

model.sol('sol6').study('std1');
model.sol('sol6').label('Parametric Solutions 2');

model.result.dataset.create('rev1', 'Revolve2D');
model.result.dataset.create('rev2', 'Revolve2D');
model.result.dataset('rev2').set('data', 'dset2');
model.result.create('pg7', 'PlotGroup2D');
model.result.create('pg6', 'PlotGroup1D');
model.result.create('pg8', 'PlotGroup3D');
model.result.create('pg9', 'PlotGroup2D');
model.result.create('pg10', 'PlotGroup3D');
model.result.create('pg11', 'PlotGroup1D');
model.result('pg7').set('data', 'dset2');
model.result('pg7').create('surf1', 'Surface');
model.result('pg7').feature('surf1').set('expr', 'meshvol');
model.result('pg6').set('data', 'dset2');
model.result('pg6').create('lngr1', 'LineGraph');
model.result('pg6').feature('lngr1').set('xdata', 'expr');
model.result('pg6').feature('lngr1').selection.set([11 12]);
model.result('pg6').feature('lngr1').set('expr', 'alpha');
model.result('pg8').create('surf1', 'Surface');
model.result('pg8').feature('surf1').set('expr', 'solid.mises');
model.result('pg8').feature('surf1').create('def', 'Deform');
model.result('pg9').create('surf1', 'Surface');
model.result('pg9').feature('surf1').set('expr', 'solid.mises');
model.result('pg9').feature('surf1').create('def', 'Deform');
model.result('pg10').create('surf1', 'Surface');
model.result('pg10').feature('surf1').set('expr', 'solid.mises');
model.result('pg10').feature('surf1').create('def', 'Deform');
model.result('pg11').create('lngr1', 'LineGraph');
model.result('pg11').feature('lngr1').set('xdata', 'expr');
model.result('pg11').feature('lngr1').set('expr', 'z');

model.study('std1').feature('param').active(false);
model.study('std1').feature('param').set('pname', {'r2' 'r1' 'dZ' 'dP'});
model.study('std1').feature('param').set('plistarr', {'500, 750, 1000' '500, 1000, 750' '6000, 4000, 3000' '100e6, 100e6, 100e6'});
model.study('std1').feature('param').set('punit', {'m' 'm' 'm' 'Pa'});

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Stationary');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('dDef').set('thresh', 0.1);
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').runAll;

model.result.dataset('rev1').set('startangle', -90);
model.result.dataset('rev1').set('revangle', 225);
model.result.dataset('rev1').set('hasspacevars', true);
model.result.dataset('rev2').set('startangle', -90);
model.result.dataset('rev2').set('revangle', 225);
model.result.dataset('rev2').set('hasspacevars', true);
model.result('pg7').label('Stress (solid)');
model.result('pg7').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg7').feature('surf1').set('colortable', 'RainbowClassic');
model.result('pg7').feature('surf1').set('smooth', 'internal');
model.result('pg7').feature('surf1').set('resolution', 'normal');
model.result('pg6').set('xlabel', 'Principal stress direction 1, z component');
model.result('pg6').set('xlabelactive', false);
model.result('pg6').feature('lngr1').set('unit', [native2unicode(hex2dec({'00' 'b0'}), 'unicode') ]);
model.result('pg6').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg6').feature('lngr1').set('xdataexpr', 'solid.sp1z');
model.result('pg6').feature('lngr1').set('xdataunit', '');
model.result('pg6').feature('lngr1').set('xdatadescr', 'Principal stress direction 1, z component');
model.result('pg6').feature('lngr1').set('smooth', 'internal');
model.result('pg6').feature('lngr1').set('resolution', 'normal');
model.result('pg8').label('Stress, 3D (solid)');
model.result('pg8').set('data', 'rev2');
model.result('pg8').set('looplevel', [3]);
model.result('pg8').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg8').feature('surf1').set('colortable', 'RainbowClassic');
model.result('pg8').feature('surf1').set('smooth', 'internal');
model.result('pg8').feature('surf1').set('resolution', 'normal');
model.result('pg8').feature('surf1').feature('def').set('revcoordsys', 'cylindrical');
model.result('pg8').feature('surf1').feature('def').set('expr', {'u' '0' 'w'});
model.result('pg8').feature('surf1').feature('def').set('descractive', true);
model.result('pg8').feature('surf1').feature('def').set('descr', 'Displacement field (material and geometry frames)');
model.result('pg9').label('Stress (solid) 1');
model.result('pg9').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg9').feature('surf1').set('colortable', 'RainbowClassic');
model.result('pg9').feature('surf1').set('smooth', 'internal');
model.result('pg9').feature('surf1').set('resolution', 'normal');
model.result('pg9').feature('surf1').feature('def').set('descr', 'Displacement field (material and geometry frames)');
model.result('pg9').feature('surf1').feature('def').set('scale', 1722.6835695782806);
model.result('pg9').feature('surf1').feature('def').set('scaleactive', false);
model.result('pg10').label('Stress, 3D (solid) 1');
model.result('pg10').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg10').feature('surf1').set('colortable', 'RainbowClassic');
model.result('pg10').feature('surf1').set('smooth', 'internal');
model.result('pg10').feature('surf1').set('resolution', 'normal');
model.result('pg10').feature('surf1').feature('def').set('revcoordsys', 'cylindrical');
model.result('pg10').feature('surf1').feature('def').set('expr', {'u' '0' 'w'});
model.result('pg10').feature('surf1').feature('def').set('descractive', true);
model.result('pg10').feature('surf1').feature('def').set('descr', 'Displacement field (material and geometry frames)');
model.result('pg11').set('xlabel', 'First principal stress (N/m<sup>2</sup>)');
model.result('pg11').set('ylabel', 'z-coordinate (m)');
model.result('pg11').set('xlabelactive', false);
model.result('pg11').set('ylabelactive', false);
model.result('pg11').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg11').feature('lngr1').set('xdataexpr', 'solid.sp1');
model.result('pg11').feature('lngr1').set('xdataunit', 'N/m^2');
model.result('pg11').feature('lngr1').set('xdatadescr', 'First principal stress');
model.result('pg11').feature('lngr1').set('smooth', 'internal');
model.result('pg11').feature('lngr1').set('resolution', 'normal');

out = model;
