clear all
file_name = 'test2dFourquad';
load_file = strcat('./tests/',file_name);
load(load_file)

FemProblem = FEM.create(file_name);
MaterialInterpolationProperties.rho_plus = 1;
MaterialInterpolationProperties.rho_minus = 0;
MaterialInterpolationProperties.E_plus = 1;
MaterialInterpolationProperties.E_minus = 1e-3;
MaterialInterpolationProperties.nu_plus = 1/3;
MaterialInterpolationProperties.nu_minus = 1/3;
MaterialInterpolation = Material_Interpolation.create(MaterialInterpolationProperties,'ISOTROPIC','SIMPALL','2D');



density = ones(FemProblem.geometry.interpolation.nelem,1);
MaterialProperties = MaterialInterpolation.computeMatProp(density);
FemProblem.preProcess();
FemProblem.setMatProps(MaterialProperties);
FemProblem.computeVariables();

if norm(FemProblem.element.K(:) - K) < 1e-12
    cprintf('green',strcat(file_name,' PASSED\n'));
else
    cprintf('err',strcat(file_name,' FAILED\n'));
end



density([1 4]) = 0.01;
MaterialProperties = MaterialInterpolation.computeMatProp(density);
FemProblem.setMatProps(MaterialProperties);
FemProblem.computeVariables();

if norm(FemProblem.element.K(:) - KwithVoidMaterial) < 1e-12
    cprintf('green',strcat(file_name,' PASSED\n'));
else
    cprintf('err',strcat(file_name,' FAILED\n'));
end