% My example
syms theta exx eyy exy real
assume(theta,'real');

% Define trigonometric terms
c2 = cos(2*theta);
s2 = sin(2*theta);

C = sym('C', [3 3]);
C(2,1) = C(1,2);
C(3,1) = 0;
C(1,3) = 0;
C(3,2) = 0;
C(2,3) = 0;
% C(1,1) = C(1,2) + 2*C(3,3);
C(2,2) = C(1,1);
assume(C, 'real');  % Optional

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Rsig matrix
Rsig = [ (1 + c2)/2, (1 - c2)/2,   s2;
         (1 - c2)/2, (1 + c2)/2,  -s2;
         -s2/2,       s2/2,       c2 ];

% Define Reps matrix
Reps = [ (1 + c2)/2, (1 - c2)/2, -s2/2;
         (1 - c2)/2, (1 + c2)/2,  s2/2;
         s2,         -s2,         c2 ];

dRsig = diff(Rsig,theta);
dReps = diff(Reps,theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_vec = [exx; eyy; exy];
rotEps = Reps*eps_vec;
rotSig = C*rotEps; 
sig = dRsig*rotSig; 
dE1 = simplify(eps_vec'*sig);

rotEps = dReps*eps_vec;
rotSig = C*rotEps; 
sig = Rsig*rotSig; 
dE2 = simplify(eps_vec'*sig);

dE = simplify(dE1 + dE2);
