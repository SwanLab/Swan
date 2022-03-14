%% Test LHS Anisotropic Stiffness

clear
clc

load('AnisotropicData.mat')
s.type = 'AnisotropicStiffnessMatrix';
s.dim.ndim = 1;
s.dim.ndimField = 1;

C = zeros(s.mesh.nnode,s.mesh.nnode,s.mesh.nelem);
for i = 1:s.mesh.nelem
    C(:,:,i) = eye(s.mesh.nnode); % Per què 4x4 en comptes de 8x8?
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