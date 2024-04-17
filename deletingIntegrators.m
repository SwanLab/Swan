%% Goal: delete integrators
% Previous goal: delete BMatrixComputer
clc; clear; close all;

% Test and trial
mesh = UnitTriangleMesh(10,10);
quad = Quadrature.set(mesh.type);
quad.computeQuadrature('QUADRATIC');
xV = quad.posgp;

trial = TrialFunction.create(mesh, 2, 'P1');
test  = TestFunction.create(mesh, 2, 'P1');

shTr = Shape(trial).evaluate(xV);
shTe = Shape(test).evaluate(xV);
dVolu(1,:,:)  = mesh.computeDvolume(quad);
mult = pagemtimes(shTr,shTe);
integ = bsxfun(@times, mult,dVolu);

M = computeElementalLHS(quad,mesh,test,trial);

shDerTe = ShapeDer(test).evaluate(xV);

function lhs = computeElementalLHS(quad,mesh, test, trial)
    xV   = quad.posgp;
    shapesTest  = test.computeShapeFunctions(xV);
    shapesTrial = trial.computeShapeFunctions(xV);
    dVolu  = mesh.computeDvolume(quad);

    nGaus  = quad.ngaus;
    nElem  = size(dVolu,2);
    nNodeTest  = size(shapesTest,1);
    nNodeTrial = size(shapesTrial,1);
    nDofTest   = nNodeTest*test.ndimf;
    nDofTrial  = nNodeTrial*trial.ndimf;

    M = zeros(nDofTest, nDofTrial, nElem);
    for igauss = 1 :nGaus
        for inode= 1:nNodeTest
            for jnode= 1:nNodeTrial
                for iunkn= 1:test.ndimf
               %     for junkn= 1:obj.trial.ndimf
                        idof = test.ndimf*(inode-1)+iunkn;
                        jdof = trial.ndimf*(jnode-1)+iunkn;
                        dvol = dVolu(igauss,:);
                        Ni = shapesTest(inode,igauss,:);
                        Nj = shapesTrial(jnode,igauss,:);
                        v = squeeze(Ni.*Nj);
                        M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                            + v(:).*dvol';
               %     end
                end
            end
        end
    end
    lhs = M;

end