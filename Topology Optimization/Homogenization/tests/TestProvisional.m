function TestProvisional
mu = sym('mu','positive');
lambda = sym('lambda','real');
direction = sym('e',[3 1],'real');
E = mu*(3*lambda+2*mu)/(lambda+mu);
nu =  lambda/(2*(lambda + mu));
A = IsotropicConstitutiveTensor3D(E,nu);



C = sym('C',[6 6],'real');
C(2,1) = C(1,2);
C(3,1) = C(1,3);
C(4,1) = C(1,4);
C(5,1) = C(1,5);
C(6,1) = C(1,6);

C(3,2) = C(2,3);
C(4,2) = C(2,4);
C(5,2) = C(2,5);
C(6,2) = C(2,6);

C(4,3) = C(3,4);
C(5,3) = C(3,5);
C(6,3) = C(3,6);

C(5,4) = C(4,5);
C(6,4) = C(4,6);

C(6,5) = C(5,6);


ft = fourthOrderTensor();
%ft.tensor = C;
ft.tensorVoigt = C;%ft.RepresentTensorinVoigt(ft.tensor);
ft.tensorVoigtInPlaneStress = ft.transform3D_2_PlaneStressInVoigt(ft.tensorVoigt);

Comp = AnisotropicContributionTensor(A,direction);


%Comp.tensorVoigt
% Comp.FirstTensor.tensorVoigt
%Comp.SecondTensor.tensorVoigt
%Comp.ThirdTensor.tensorVoigt


Mod2D = simplify(Comp.tensorVoigtInPlaneStress);
Iso2D = simplify(Comp.FirstTensor.tensorVoigtInPlaneStress);
Ani2D1 = Comp.SecondTensor.tensorVoigtInPlaneStress;
ANi2D2 = Comp.ThirdTensor.tensorVoigtInPlaneStress;


lambda = simplify(Comp.FirstTensor.lambda);
mu = simplify(Comp.FirstTensor.lambda);

PSnotation = subs(Comp.tensorVoigtInPlaneStress,'lambda',2*mu*lambda/(2*mu-lambda));
New = simplify(PSnotation)
end