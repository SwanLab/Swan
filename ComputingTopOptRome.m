clear
clc
close all

% INPUTS
x1                       = 2;
y1                       = 1;
N                        = 200; % 600
M                        = 100; % 300
s.testName               = 'romeCase.m';
s.filename               = "'archRome'";
s.ptype                  = "'MACRO'";
s.initial_case           = "'full'";
s.cost                   = "'compliance'";
s.weights                = '[1]'; % 0.25
s.optimizerUnconstrained = "'SLERP'";
s.designVariable         = "'LevelSet'";
s.nsteps                 = '1';
s.Vfrac_final            = '0.4';
s.Vfrac_initial          = '1';

% END INPUTS

testWriter(s);

filename      = char(s.filename);
filename      = filename(2:end-1);
filename      = [char(filename),'.m'];
z.testName    = filename;
z.problemCase = filename(1:end-6);
z.x1          = x1;
z.y1          = y1;
z.N           = N;
z.M           = M;
z.P           = -1;
z.DoF         = 2;
Benchmark = FEMInputWriter(z);
Benchmark.createTest();

ss.testName = s.testName;
t = TopOptComputer(ss);
k.aJmax = 2;
k.aGmax = 0.01;
t.compute(k);

% Pending:
% - Find Dehomogenizer - saveimage - obtain gif single simulation
% - Big loop over properties and cases, and obtain gif for different simus compliance (then obtain aJ,aG callibration)
% - Jump to compliance + different kind of perimeters (20 steps) [use epsilon y exp por defecto]


function testWriter(cParams)

fileID = fopen(['Input/',cParams.testName],'w');
fprintf(fileID,['filename = ',char(cParams.filename),';\n']);
fprintf(fileID,['ptype = ',char(cParams.ptype),';\n']);
fprintf(fileID,['method = ',char("'SIMPALL'"),';\n']);
fprintf(fileID,['materialType = ',char("'ISOTROPIC'"),';\n']);
fprintf(fileID,['initial_case = ',char(cParams.initial_case),';\n']);
fprintf(fileID,['cost = {',char(cParams.cost),'};\n']);
fprintf(fileID,['weights =',cParams.weights,';\n']);
fprintf(fileID,['constraint = {',char("'volumeConstraint'"),'};\n']);
fprintf(fileID,['constraint_case = {',char("'EQUALITY'"),'};\n']);
fprintf(fileID,['optimizerUnconstrained = ',char(cParams.optimizerUnconstrained),';\n']);
fprintf(fileID,['optimizer = ',char("'NullSpace'"),';\n']);
fprintf(fileID,['incrementFactor = 2',';\n']);
fprintf(fileID,['designVariable = ',char(cParams.designVariable),';\n']);
fprintf(fileID,['filterType = ',char("'P1'"),';\n']);
fprintf(fileID,['nsteps = ',cParams.nsteps,';\n']);
fprintf(fileID,['Vfrac_final = ',cParams.Vfrac_final,';\n']);
fprintf(fileID,'optimality_final = 1e-3;\n');
fprintf(fileID,'constr_final =1e-3;\n');
fprintf(fileID,['Vfrac_initial = ',cParams.Vfrac_initial,';\n']);
fprintf(fileID,'optimality_initial = 1e-3;\n');
fprintf(fileID,'constr_initial = 1e-3;\n');
fprintf(fileID,'Perimeter_target = 5;\n');
fprintf(fileID,'perimeterTarget = 5;\n');
fprintf(fileID,'TOL.rho_plus = 1;\n');
fprintf(fileID,'TOL.rho_minus = 0;\n');
fprintf(fileID,'TOL.E_plus = 1;\n');
fprintf(fileID,'TOL.E_minus = 1e-3;\n');
fprintf(fileID,'TOL.nu_plus = 1/3;\n');
fprintf(fileID,'TOL.nu_minus = 1/3;\n');
fprintf(fileID,'plotting = true;\n');
fprintf(fileID,'printing = false;\n');
fprintf(fileID,'printing_physics = false;\n');
fprintf(fileID,'monitoring = true;\n');
fprintf(fileID,'monitoring_interval = 1;\n');
fprintf(fileID,'maxiter = 300*nsteps;');
fclose(fileID);

end