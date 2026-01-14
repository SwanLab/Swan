function wNEW = WeightsHomogenization(BasisPone,BasisPoneZ,ngaus,nstrain,W,wINTF,VOL) 

if nargin == 0
    load('tmp1.mat')
end

% Construct matrix Y 
% ------------------

Y = MatrixBasisStressesReshape(BasisPone,nstrain,length(W)) ; 
%[UU,SS,VV] = RSVDT(Y) ; 

% Exact Integral 
b = Y'*W ; 
% Integral Using the weights computed for internal forces
Yz = MatrixBasisStressesReshape(BasisPoneZ,nstrain,length(wINTF)) ; 

bWintF = Yz'*wINTF ; 

errorWeightsINTF= norm(b-bWintF)/norm(b) ; 

disp('**********************************************************************************')
disp(['ERROR homog. using intern. force. weights = ',num2str(errorWeightsINTF)]) ; 
disp('**********************************************************************************')
wNEW = lsqnonneg(Yz',b) ; 
%[UU,SS,VV] = RSVDT(Yz) ;

bTAILORED = Yz'*wNEW ; 
errorNEW= norm(b-bTAILORED)/norm(b) ; 
disp('**********************************************************************************')
disp(['ERROR homog. using tailored  weights = ',num2str(errorNEW)]) ; 
disp('**********************************************************************************')

wNEW = wNEW/VOL ; 

end

function Y = MatrixBasisStressesReshape(BasisPone,nstrain,ngausTOTAL)


ncols = size(BasisPone,2)*nstrain ; 
Y = zeros(ngausTOTAL,ncols) ; 
iniCOLnew = 1; 
for istrain = 1:nstrain 
    ROWS_old = istrain:nstrain:size(BasisPone,1) ; 
    finROWnew = iniCOLnew + size(BasisPone,2)-1 ;      
    Y(:,iniCOLnew:finROWnew)    =  BasisPone(ROWS_old,:) ; 
    iniCOLnew = finROWnew + 1; 
end


end
 