%% Test LHS Anisotropic Stiffness

clear
clc

load('AnisotropicData.mat')
s.type = 'AnisotropicStiffnessMatrix';
s.dim.ndim = 1;
s.dim.ndimField = 1;
s.dofsInElem = s.globalConnec';
% s.material = 'anisotropy';

C = zeros(s.mesh.ndim,s.mesh.ndim,s.mesh.nelem);
for i = 1:s.mesh.nelem
    C(:,:,i) = eye(s.mesh.ndim);
end
s.Celas = C;

LHSint = LHSintegrator.create(s);
LHStest = LHSint.compute();

error = full(max(max(abs(LHS-LHStest))));

if error<1e-6
    disp('TestPassed')
else
    disp('TestFailed')
end