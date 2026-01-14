function stressEL = Stress_Elements_Average(weig,stress,ngausE,nstress,nelem)
% ------------------------------------------------
%%% Average stresses ON ELEMENTS 
% ------------------------------------------------
%weig = DATAOUT.wSTs ; % Weights at all Gauss points 
% volume of each element 
volE = reshape(weig,ngausE,[]) ; 
volE = sum(volE) ; 
% ------------------------------------------
% Reshape stress array 
stress = reshape( stress,nstress,[]) ; 
% Multiply the stress at each Gauss point by its corresponding weight 
stress = (bsxfun(@times,stress',weig))' ; 
% Calculate the sum of weighted stresses at each element 
% ----------------------------------------------
 stressEL = zeros(nstress,nelem) ; 
for igaus = 1:ngausE
    stressEL = stressEL + stress(:,igaus:ngausE:end) ; 
end
% Divide by the volume of each element 
stressEL = (bsxfun(@rdivide,stressEL',volE'))' ;  
%stressEL = stressEL(:) ; 