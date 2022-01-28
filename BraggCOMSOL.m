clear all
close all
clc
tic
%% Variabili per controllare le dimensioni della geometria
nH = 2.1306; % Tantala
nL = 1.4585; % Silica
n_cavity = 1.41;
n_air = 1;
lambda_thickness = 910; %nm
LH = lambda_thickness/(4*nH);
LL = lambda_thickness/(4*nL);
L_cavity = lambda_thickness/(2*n_cavity);
ray_Bragg = 1000; % Ray in the structure
rHSQ = 500;
h_PML = 1000;
%% Bragg Adrien
% Reference wavelength: 910nm
% 
% 1	-3	0.97086
% 2	-2	1.02367
% 3	-3	1.09909
% 4	-2	1.07653
% 5	-3	1.03703
% 6	-2	0.97555
% 7	-3	0.97214
% 8	-2	0.97688
% 9	-3	0.96318
% 10	-2	0.9915
% 11	-3	1.00896
% 12	-2	1.03276
% 13	-3	1.02587
% 14	-2	0.99642
% 15	-3	0.97038
% 16	-2	0.97342
% 17	-3	0.99916
% 18	-2	1
% 19	-3	1
% 20	-2	1
% 21	-3	1
% 22	-2	1
% 23	-3	1
% 24	-2	1
% 25	-3	1
% 26	-2	1
% 27	-3	1
% 28	-2	1
% 29	-3	1
% 30	-2	1
% 
% First column: layer number, starts at the substrate, "30" is
% topmost Second column: "-3" denotes Tantala, "-2" is Silica Third

thicknesses = [0.970860000000000;1.02367000000000;1.09909000000000;1.07653000000000;1.03703000000000;...
    0.975550000000000;0.972140000000000;0.976880000000000;0.963180000000000;0.991500000000000;...
    1.00896000000000;1.03276000000000;1.02587000000000;0.996420000000000;0.970380000000000;...
    0.973420000000000;0.999160000000000;1;1;1;1;1;1;1;1;1;1;1;1;1];

materials = [-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2;-3;-2];

for i = 1 : size(materials,1)
    
    if materials(i) == -3
        materials(i) = nH;
        thicknesses(i) = thicknesses(i)*LH;

        
    elseif materials(i) == -2
        materials(i) = nL;
        thicknesses(i) = thicknesses(i)*LL;
    end
end
thicknesses = thicknesses';
thicknesses = flip(thicknesses,2);
materials = materials';
materials = flip(materials,2);
%% Import COMSOL class
nom=['/home/rapisarda/Documents/PhD/COMSOL/','prova_bojh.mph'];
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
model.modelPath('/home/rapisarda/Documents/PhD/COMSOL');
model.modelNode.create('comp1');
%% Parametric sweep on the HSQ thickness
h_corrected = [3.95996518711923e-07,3.68227054920752e-07,3.53244059467294e-07,...
    3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,...
    3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,...
    3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,...
    3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07,3.24166429182103e-07];

radius = 300:100:2499;
eigenfreq = zeros(5,size(h_corrected,2));
Q_factor = zeros(5,size(h_corrected,2));


%% Create geometry sequence 
model.geom.create('Bragg', 2);
model.geom('Bragg').axisymmetric(true);
%% Mesh
model.mesh.create('mesh1', 'Bragg');
%% Add physics
model.physics.create('emw', 'ElectromagneticWaves', 'Bragg');
model.physics('emw').feature('wee1').set('DisplacementFieldModel', 'RefractiveIndex');
model.physics('emw').prop('outofplanewavenumber').set('mFloquet', '1');
model.physics('emw').feature('wee1').setIndex('materialType', 'solid', 0);
%% Create study
model.study.create('std1');
model.study('std1').create('eig', 'Eigenfrequency');
model.study('std1').feature('eig').activate('emw', true);
model.study('std1').feature('eig').set('shift', '2.6e14'); %%%% Per scegliere la frequenza attorno a cui cercare
model.study('std1').feature('eig').set('neigs', 5); %%% Per scegliere quante aufofrequenze cercare
% for r = 1:size(h_corrected,2)
% model.param.set('hCavity',num2str(h_corrected(r)));
% model.param.set('rHSQ',num2str(radius(r)))
%% Add global parameters 
model.param.set('rHSQ','2499');
model.param.set('lambda','910');
model.param.set('n_cavity','1.41');
model.param.set('n_air','1');
model.param.set('nH','2.1306');
model.param.set('nL','1.4585');
model.param.set('ray_Bragg','3500');
model.param.set('hCavity','lambda/(2*n_cavity)');
model.param.set('h_PML','1000');



%% Tantala sequence
model.geom('Bragg').scaleUnitValue(true);
model.geom('Bragg').lengthUnit('nm');
tag = strcat('Tantala',num2str(1));
Tantala = model.geom('Bragg').feature.create(tag, 'Rectangle');
Tantala.set('size', [mphevaluate(model,'ray_Bragg'),thicknesses(1)]);
Tantala.set('pos',[0 0]);

for i = 3:2:size(materials,2)
tag = strcat('Tantala',num2str(i));
Tantala = model.geom('Bragg').feature.create(tag, 'Rectangle');
Tantala.set('size', [mphevaluate(model,'ray_Bragg'),thicknesses(i)]);
Tantala.set('pos',[0 sum(thicknesses(1:i-1))]);
end
%% Silica sequence
tag = strcat('Silica',num2str(2));
Silica = model.geom('Bragg').feature.create(tag, 'Rectangle');
Silica.set('size', [mphevaluate(model,'ray_Bragg'),thicknesses(2)]);
Silica.set('pos',[0 thicknesses(1)]);

for j = 4:2:size(materials,2)
tag = strcat('Silica',num2str(j));
Silica = model.geom('Bragg').feature.create(tag, 'Rectangle');
Silica.set('size', [mphevaluate(model,'ray_Bragg'),thicknesses(j)]);
Silica.set('pos',[0 sum(thicknesses(1:j-1))]);
end
%% Cavity sequence
tag = 'Air';
Air = model.geom('Bragg').feature.create(tag, 'Rectangle');
Air.set('size', [mphevaluate(model,'ray_Bragg'),mphevaluate(model,'hCavity')]);
Air.set('pos',[0 mphevaluate(model,'hCavity')*(-1)]);

tag = 'HSQ';
HSQ = model.geom('Bragg').feature.create(tag, 'Rectangle');
HSQ.set('size', {'rHSQ' 'hCavity'});
HSQ.set('pos',[0 mphevaluate(model,'hCavity')*(-1)]);

%% PML sequence
tag = strcat('PML_SIDE',num2str(2));
PML_SIDE = model.geom('Bragg').feature.create(tag, 'Rectangle');
PML_SIDE.set('size', [1000,mphevaluate(model,'hCavity')+2*sum(thicknesses)]);
PML_SIDE.set('pos',[mphevaluate(model,'ray_Bragg')-1000 (sum(thicknesses))*(-1)-mphevaluate(model,'hCavity')]);

tag = strcat('PML_TOP',num2str(2));
PML_TOP = model.geom('Bragg').feature.create(tag, 'Rectangle');
PML_TOP.set('size', [mphevaluate(model,'ray_Bragg'),1000]);
PML_TOP.set('pos',[0 sum(thicknesses)]);

tag = strcat('PML_BOTTOM',num2str(2));
PML_BOTTOM = model.geom('Bragg').feature.create(tag, 'Rectangle');
PML_BOTTOM.set('size', [mphevaluate(model,'ray_Bragg'),1000]);
PML_BOTTOM.set('pos',[0 sum(thicknesses)*(-1)+mphevaluate(model,'hCavity')*(-1)-mphevaluate(model,'h_PML')]);

%% Bottom Bragg

model.geom('Bragg').create('copy1', 'Copy');
model.geom('Bragg').feature('copy1').selection('input').set({'Tantala1' 'Tantala3' 'Tantala5' 'Tantala7' 'Tantala9'...
    'Tantala11' 'Tantala13' 'Tantala15' 'Tantala17' 'Tantala19' 'Tantala21' 'Tantala23' 'Tantala25'...
    'Tantala27' 'Tantala29' 'Silica2' 'Silica4' 'Silica6' 'Silica8' 'Silica10' 'Silica12' 'Silica14' 'Silica16'...
    'Silica18' 'Silica20' 'Silica22' 'Silica24' 'Silica26' 'Silica28' 'Silica30'});
model.geom('Bragg').feature('copy1').set('disply', sum(thicknesses)*(-1)+mphevaluate(model,'hCavity')*(-1) );


%% Rotate the Bottom Bragg to have a symmetric cavity
model.geom('Bragg').create('rot1', 'Rotate');
model.geom('Bragg').feature('rot1').selection('input').set({'copy1'});
model.geom('Bragg').feature('rot1').set('rot', '180');
model.geom('Bragg').feature('rot1').set('pos', {num2str((mphevaluate(model,'ray_Bragg'))/2) num2str(sum(thicknesses)/2*(-1)+mphevaluate(model,'hCavity')*(-1))});

model.geom('Bragg').run;
%% Insert PML
model.coordSystem.create('pml1', 'Bragg', 'PML');
model.coordSystem('pml1').selection.set([1 63 65:125]); 
model.coordSystem('pml1').set('ScalingType', 'Cylindrical');
%% Create Materials
% Air
model.material.create('mat1', 'Common', 'comp1');
model.material('mat1').label('Air');
model.material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat1').propertyGroup('RefractiveIndex').set('n', mphevaluate(model,'n_air'));
model.material('mat1').selection.set([1 63 64 95]);

% HSQ
model.material.create('mat2', 'Common', 'comp1');
model.material('mat2').label('HSQ');
model.material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat2').propertyGroup('RefractiveIndex').set('n', mphevaluate(model,'n_cavity'));
model.material('mat2').selection.set(32);

% Silica
model.material.create('mat3', 'Common', 'comp1');
model.material('mat3').label('Silica');
model.material('mat3').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat3').propertyGroup('RefractiveIndex').set('n', mphevaluate(model,'nL'));
model.material('mat3').selection.set([3:2:61 66:2:124]);

% Tantala
model.material.create('mat4', 'Common', 'comp1');
model.material('mat4').label('Tantala');
model.material('mat4').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat4').propertyGroup('RefractiveIndex').set('n', mphevaluate(model,'nH'));
model.material('mat4').selection.set([2:2:30 34:2:62 65:2:93 97:2:125 ]);
%% To show on MatLAB the geometry that we built
mphgeom(model,'Bragg');

%% Studio
model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('eig').set('notlistsolnum', 1);
model.study('std1').feature('eig').set('notsolnum', '1');
model.study('std1').feature('eig').set('listsolnum', 1);
model.study('std1').feature('eig').set('solnum', '1');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'eig');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'eig');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('control', 'eig');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').attach('std1');

model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').label('Electric Field (emw)');
model.result('pg1').set('oldanalysistype', 'noneavailable');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').set('data', 'dset1');
model.result('pg1').feature.create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('oldanalysistype', 'noneavailable');
model.result('pg1').feature('surf1').set('data', 'parent');
model.result.dataset.create('rev1', 'Revolve2D');
model.result.dataset('rev1').set('startangle', -90);
model.result.dataset('rev1').set('revangle', 225);
model.result.dataset('rev1').set('data', 'dset1');

model.sol('sol1').runAll;

model.result('pg1').run;

model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('gev1').set('expr', 'emw.Qfactor');
model.result.numerical('gev1').set('descr', 'Quality factor');
model.result.table.create('tbl1', 'Table');
model.result.table('tbl1').comments('Global Evaluation 1 (emw.Qfactor)');
model.result.numerical('gev1').set('table', 'tbl1');
model.result.numerical('gev1').setResult;

mphsave(model,nom)

%% Display data
eval = mpheval(model,'emw.freq','Dataset','dset1','edim',0,'selection',1);
eigenfreq(:,1) = eval.d1; 

Q = mpheval(model,'emw.Qfactor','Dataset','dset1','edim',0,'selection',1);
Q_factor(:,1) = Q.d1;

for i = 1:size(eigenfreq,1)
figure(i)
stringa = strcat('with(',num2str(i),',emw.normE)');
model.result('pg1').feature('surf1').set('expr', stringa);
mphplot(model, 'pg1','rangenum',1);
hold on 
mphgeom(model,'Bragg','facemode','off');
hold off
end
% stampa = strcat('Raggio appena analizzato:',num2str(radius(r)));
% disp(stampa);
% end
toc
% tt = load ('handel');
% sound(tt.y)