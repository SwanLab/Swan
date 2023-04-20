clear
clc
close all

x1vec = [2;6;2];
y1vec = [1;1;1];
Nvec  = 0.5*[200;600;200];
Mvec  = 0.5*[100;100;100];
filenamevec = ["'archRome'";"'bridgeRome'";"'cantileverRome'"];
VfracFinalVec = ['0.15';'0.50';'0.50'];
aJvec =[4.5;5;3.5];
aGvec = [0.15;0.02;0.06];

costVec = ["'compliance','anisotropicPerimeter2D'";"'compliance','anisotropicPerimeterInterior2D'"]; % cambiar a mano el nu y alpha

optUnconstrVec = ["'SLERP'";"'PROJECTED GRADIENT'"];
designVarVec   = ["'LevelSet'";"'Density'"];


for i=2:2%length(x1vec)
    for j=2:2%size(costVec,1)
        for l=2:2%size(optUnconstrVec,1)
            close all
            % INPUTS
            x1                       = x1vec(i);
            y1                       = y1vec(i);
            N                        = Nvec(i); % 600 arch
            M                        = Mvec(i); % 300 arch
            s.testName               = ['romeCase',int2str(i),int2str(j),int2str(l)];
            s.filename               = filenamevec(i);
            s.ptype                  = "'MACRO'";
            s.initial_case           = "'full'";
            s.cost                   = costVec(j);
            s.weights                = '[1,0.18]'; % 0.25
            s.optimizerUnconstrained = optUnconstrVec(l);
            s.designVariable         = designVarVec(l);
            s.nsteps                 = '10';
            s.Vfrac_final            = VfracFinalVec(i,:);
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
            z.P           = -10;
            z.DoF         = 2;
            Benchmark = FEMInputWriter(z);
            Benchmark.createTest();

            ss.testName = s.testName;
            t = TopOptComputer(ss);
            k.aJmax = aJvec(i);
            k.aGmax = aGvec(i);
            t.compute(k);

            % - Find Dehomogenizer - saveimage - obtain gif simulations LS and density

        end
    end
end

function testWriter(cParams)

fileID = fopen(['Input/',cParams.testName,'.m'],'w');
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
fprintf(fileID,'plotting = false;\n');
fprintf(fileID,'printing = true;\n');
fprintf(fileID,'printing_physics = false;\n');
fprintf(fileID,'monitoring = true;\n');
fprintf(fileID,'monitoring_interval = 1;\n');
fprintf(fileID,'maxiter = 200*nsteps;');
fclose(fileID);

end