function [aFINAL,NAMES ]= DesignTestPlates(DISP,ROT,Ly,Lx)
% 

if nargin == 0
    DISP = 0.1; 
    ROT = 0.1 ; 
    Ly = 1 ; 
    Lx  = 1; 
end


NAMES{1} = 'StrainXX' ;  %Trans x, face 3
NAMES{2} = 'StrainXY' ;  % Trans y, face 3
NAMES{3} = 'StrainXZ' ;  % Trans z, face 3
NAMES{4} = 'TorsionX' ;  % Rot x, face 3
NAMES{5} = 'BendingXZ' ;  % Rot y, face 3
NAMES{6} = 'BendingXY' ;  % Rot z, face 3
%%% 
NAMES{7} = 'StrainY' ;  % Trans y, face 4
NAMES{8} = 'StrainYZ' ; % Trans z, face 4
NAMES{9} = 'BendingYZ' ; % Rot X, face 4
NAMES{10} = 'TorsionY' ; % Rot Y, face 4 
NAMES{11} = 'BendingYX' ; % Roz Z, face 4 
%%%%%  
NAMES{12} = 'Null1' ;  % Rot x, face 1
% 
NAMES{13} = 'Null2'; 
NAMES{14}  ='Null3' ; 
 


% RIGID BODY MODES 
% ----------------
aRB = RigidBody(Lx,Ly,DISP,ROT) ; 

%% TESTS moving FACE 3 
aFACE3 = Face3(Lx,Ly,DISP,ROT) ; 
%  TESTS moving FACE 4
aFACE4 = Face4(Lx,Ly,DISP,ROT) ; 

% Check whether motions aFACE4 are linearly depedent from aRB and aFACE3 
aRB3 = [aRB,aFACE3] ; 
aTRAJ = [aFACE3]; 
for i=1:size(aFACE4,2)
   rB3_ad =   rank([aRB3,aFACE4(:,i)]) ;
   if rB3_ad == (size(aRB3,2)+1)
        aTRAJ = [aTRAJ,aFACE4(:,i)] ; 
   end
end

 

% %  TESTS moving FACE 2
% aRB3 = [aTRAJ,aRB] ; 
% aFACE2 = Face2(Lx,Ly,DISP,ROT) ; 
% 
% % Adding new modes
% for i=1:size(aFACE2,2)
%    rB3_ad =   rank([aRB3,aFACE2(:,i)]) ;
%    if rB3_ad == (size(aRB3,2)+1)
%         aTRAJ = [aTRAJ,aFACE2(:,i)] ; 
%    end
% end






%%% Now we seek the orthogonal complement 
MMM = max(max(aTRAJ)); 
aCOMP = null([aTRAJ,aRB]'); 

for inull = 1:size(aCOMP,2)
    NNN = max(max(aCOMP(:,inull))); 
    aCOMP(:,inull) = aCOMP(:,inull)*MMM/NNN ;  
end

aFINAL = [aTRAJ,aCOMP] ; 




end

function  aFACE3 = Face3(Lx,Ly,DISP,ROT)
 

%% TESTS moving FACE 3 
% -------------
nmodes = 6 ; 
nDOF  = 5; 
a{1} = zeros(nDOF,nmodes) ; 
a{2} = zeros(nDOF,nmodes) ; 
a{3} = zeros(nDOF,nmodes) ; 
a{4} = zeros(nDOF,nmodes) ; 

 
% -------------------------------
% Test 1  (normal strain x)
% 
itest = 1; 
a{3}(1,itest) = DISP ; 
a{4}(1,itest) = DISP ; 
% Test 2   (in-plane shear)
itest  = 2; 
a{3}(2,itest) = DISP ; 
a{4}(2,itest) = DISP ; 
% Test 3   (out-of-plane shear)
itest  = 3; 
a{3}(3,itest) = DISP ; 
a{4}(3,itest) = DISP ; 
% Test 4   (torsion)
itest  = 4; 
a{3}(4,itest) = ROT ; 
a{3}(3,itest) = ROT*Ly/2 ; 
a{4}(4,itest) = ROT ; 
a{4}(3,itest) = -ROT*Ly/2 ; 
% Test 5   (out-of-plane bending)
itest  = 5; 
a{3}(5,itest) = ROT ; 
a{4}(5,itest) = ROT ;
% Test 6 (in-plane bending)
itest  = 6; 
a{3}(1,itest) = ROT*Ly/2 ; 
a{4}(1,itest) = -ROT*Ly/2 ; 


aFACE3 = cell2mat(a') ; 


end
 
 
function  aFACE4 = Face4(Lx,Ly,DISP,ROT)
 

%% TESTS moving FACE 4
% -------------
nmodes = 6 ; 
nDOF  = 5; 
a{1} = zeros(nDOF,nmodes) ; 
a{2} = zeros(nDOF,nmodes) ; 
a{3} = zeros(nDOF,nmodes) ; 
a{4} = zeros(nDOF,nmodes) ; 

 
% -------------------------------
% Test 1  (in-plane shear )
% 
itest = 1; 
a{1}(1,itest) = DISP ; 
a{4}(1,itest) = DISP ; 
% Test 2   (normal strain)
itest  = 2; 
a{1}(2,itest) = DISP ; 
a{4}(2,itest) = DISP ; 
% Test 3   (out-of-plane shear)
itest  = 3; 
a{1}(3,itest) = DISP ; 
a{4}(3,itest) = DISP ; 
% Test 5   (torsion)
itest  = 5; 
a{1}(5,itest) = ROT ; 
a{1}(3,itest) = -ROT*Lx/2 ; 
a{4}(5,itest) = ROT ; 
a{4}(3,itest) = +ROT*Lx/2 ; 
% Test 4   (out-of-plane bending)
itest  = 4; 
a{1}(4,itest) = ROT ; 
a{4}(4,itest) = ROT ;
% Test 6 (in-plane bending)
itest  = 6; 
a{1}(2,itest) = -ROT*Lx/2 ; 
a{4}(2,itest) = +ROT*Lx/2 ; 


aFACE4 = cell2mat(a') ; 


end
  

 
 

function  aRB = RigidBody(Lx,Ly,DISP,ROT)
 

%% TESTS moving FACE 4
% -------------
nmodes = 6 ; 
nDOF  = 5; 
a{1} = zeros(nDOF,nmodes) ; 
a{2} = zeros(nDOF,nmodes) ; 
a{3} = zeros(nDOF,nmodes) ; 
a{4} = zeros(nDOF,nmodes) ; 

 
% -------------------------------
% Test 1   - Translation x 
% 
itest = 1; 
a{1}(1,itest) = DISP ; 
a{2}(1,itest) = DISP ; 
a{3}(1,itest) = DISP ; 
a{4}(1,itest) = DISP ; 
% Test 2     - Translation y 
itest = 2 ; 
a{1}(2,itest) = DISP ; 
a{2}(2,itest) = DISP ; 
a{3}(2,itest) = DISP ; 
a{4}(2,itest) = DISP ; 
% Test 3     - Translation z
itest = 3 ; 
a{1}(3,itest) = DISP ; 
a{2}(3,itest) = DISP ; 
a{3}(3,itest) = DISP ; 
a{4}(3,itest) = DISP ; 
% Test 4
itest = 4 ; 
a{1}(4,itest) = ROT ; 
a{2}(4,itest) = ROT ; 
a{3}(4,itest) = ROT ; 
a{4}(4,itest) = ROT ; 

a{1}(3,itest) = ROT*Lx/2 ;    
a{2}(3,itest) = -ROT*Lx/2 ;  
a{3}(3,itest) = -ROT*Lx/2  ;  
a{4}(3,itest) = ROT*Lx/2 ;   

% Test 5
itest = 5 ; 
a{1}(5,itest) = ROT ; 
a{2}(5,itest) = ROT ; 
a{3}(5,itest) = ROT ; 
a{4}(5,itest) = ROT ; 

a{1}(3,itest) = ROT*Ly/2 ;    
a{2}(3,itest) = ROT*Ly/2 ;  
a{3}(3,itest) = -ROT*Ly/2  ;  
a{4}(3,itest) = -ROT*Ly/2 ;   


% Test 6
itest = 6 ; 
a{1}(1,itest) = -ROT*Ly/2 ;    
a{2}(1,itest) = ROT*Ly/2 ;  
a{3}(1,itest) = ROT*Ly/2  ;  
a{4}(1,itest) = -ROT*Ly/2 ;  

a{1}(2,itest) = -ROT*Lx/2 ;    
a{2}(2,itest) = -ROT*Lx/2 ;  
a{3}(2,itest) =+ROT*Lx/2  ;  
a{4}(2,itest) = +ROT*Lx/2 ;  

aRB = cell2mat(a') ; 


end
 







