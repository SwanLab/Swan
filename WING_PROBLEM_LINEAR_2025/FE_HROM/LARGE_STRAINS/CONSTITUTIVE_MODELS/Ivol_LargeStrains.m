function Ivol = Ivol_LargeStrains(Cb,nstrain) ;

if nargin == 0
    load('tmp.mat')
end

Ivol = zeros(size(Cb,1),nstrain) ;


SROWS = cell(1,nstrain) ;
for irows  =1:nstrain
    SROWS{irows} = irows:nstrain:size(Cb,1) ;
end

if nstrain == 3
    % Generated automatically by
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/CONSTITUTIVE_MODELS/Ivol_LargeStrains_aux.m
 Ivol(SROWS{1},1) = Cb(SROWS{1}).^2; 
Ivol(SROWS{1},2) = Cb(SROWS{1}).*Cb(SROWS{2}); 
Ivol(SROWS{1},3) = Cb(SROWS{1}).*Cb(SROWS{3}); 
Ivol(SROWS{2},1) = Ivol(SROWS{1},2); 
Ivol(SROWS{2},2) = Cb(SROWS{2}).^2; 
Ivol(SROWS{2},3) = Cb(SROWS{2}).*Cb(SROWS{3}); 
Ivol(SROWS{3},1) = Ivol(SROWS{1},3); 
Ivol(SROWS{3},2) = Ivol(SROWS{2},3); 
Ivol(SROWS{3},3) = Cb(SROWS{3}).^2; 

    
elseif nstrain == 6
    % Generated automatically by
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/CONSTITUTIVE_MODELS/Ivol_LargeStrains_aux.m
   Ivol(SROWS{1},1) = Cb(SROWS{1}).^2; 
Ivol(SROWS{1},2) = Cb(SROWS{1}).*Cb(SROWS{2}); 
Ivol(SROWS{1},3) = Cb(SROWS{1}).*Cb(SROWS{3}); 
Ivol(SROWS{1},4) = Cb(SROWS{1}).*Cb(SROWS{4}); 
Ivol(SROWS{1},5) = Cb(SROWS{1}).*Cb(SROWS{5}); 
Ivol(SROWS{1},6) = Cb(SROWS{1}).*Cb(SROWS{6}); 
Ivol(SROWS{2},1) = Ivol(SROWS{1},2); 
Ivol(SROWS{2},2) = Cb(SROWS{2}).^2; 
Ivol(SROWS{2},3) = Cb(SROWS{2}).*Cb(SROWS{3}); 
Ivol(SROWS{2},4) = Cb(SROWS{2}).*Cb(SROWS{4}); 
Ivol(SROWS{2},5) = Cb(SROWS{2}).*Cb(SROWS{5}); 
Ivol(SROWS{2},6) = Cb(SROWS{2}).*Cb(SROWS{6}); 
Ivol(SROWS{3},1) = Ivol(SROWS{1},3); 
Ivol(SROWS{3},2) = Ivol(SROWS{2},3); 
Ivol(SROWS{3},3) = Cb(SROWS{3}).^2; 
Ivol(SROWS{3},4) = Cb(SROWS{3}).*Cb(SROWS{4}); 
Ivol(SROWS{3},5) = Cb(SROWS{3}).*Cb(SROWS{5}); 
Ivol(SROWS{3},6) = Cb(SROWS{3}).*Cb(SROWS{6}); 
Ivol(SROWS{4},1) = Ivol(SROWS{1},4); 
Ivol(SROWS{4},2) = Ivol(SROWS{2},4); 
Ivol(SROWS{4},3) = Ivol(SROWS{3},4); 
Ivol(SROWS{4},4) = Cb(SROWS{4}).^2; 
Ivol(SROWS{4},5) = Cb(SROWS{4}).*Cb(SROWS{5}); 
Ivol(SROWS{4},6) = Cb(SROWS{4}).*Cb(SROWS{6}); 
Ivol(SROWS{5},1) = Ivol(SROWS{1},5); 
Ivol(SROWS{5},2) = Ivol(SROWS{2},5); 
Ivol(SROWS{5},3) = Ivol(SROWS{3},5); 
Ivol(SROWS{5},4) = Ivol(SROWS{4},5); 
Ivol(SROWS{5},5) = Cb(SROWS{5}).^2; 
Ivol(SROWS{5},6) = Cb(SROWS{5}).*Cb(SROWS{6}); 
Ivol(SROWS{6},1) = Ivol(SROWS{1},6); 
Ivol(SROWS{6},2) = Ivol(SROWS{2},6); 
Ivol(SROWS{6},3) = Ivol(SROWS{3},6); 
Ivol(SROWS{6},4) = Ivol(SROWS{4},6); 
Ivol(SROWS{6},5) = Ivol(SROWS{5},6); 
Ivol(SROWS{6},6) = Cb(SROWS{6}).^2; 

    
    
else
    error('Option not implemented')
end
