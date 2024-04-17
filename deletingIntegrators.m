%% Goal: delete integrators
% Previous goal: delete BMatrixComputer
clc; clear; close all;

% Test and trial
mesh = UnitTriangleMesh(10,10);
quad = Quadrature.set(mesh.type);
quad.computeQuadrature('QUADRATIC');
xV = quad.posgp;

ndimf = 2;
trial = TrialFunction.create(mesh, ndimf, 'P1');
test  = TestFunction.create(mesh, ndimf, 'P1');

shTr = Shape(trial).evaluate(xV);
shTe = Shape(test).evaluate(xV);
dVolu(1,:,:)  = mesh.computeDvolume(quad);
mult = pagemtimes(shTr,shTe');
integ = bsxfun(@times, mult, dVolu);
% kron(squeeze(A(:,:,1)),eye(2)) % !!!!
integ = openprod(integ,eye(ndimf));

M = massmatrix(quad,mesh,test,trial);

shDerTe = ShapeDer(test).evaluate(xV);

function z = openprod(A,B)
    z = bsxfun(@times, permute(A, [4 1 5 2 3]), permute(B, [1 4 2 5 3]));  %// step 1
    z = reshape(z, size(A,1)*size(B,1), size(A,2)*size(B,2), size(A,3)); 
end

function lhs = massmatrix(quad,mesh, test, trial)
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