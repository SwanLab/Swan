function indPLASm = INDEXcomponents(indPLAS,nstrain) 

indPLASm = zeros(nstrain,length(indPLAS)) ; 
for istrain=1:nstrain 
    indPLASm(istrain,:) = (indPLAS-1)*nstrain+istrain ; 
end
indPLASm = reshape(indPLASm,length(indPLAS)*nstrain,1);