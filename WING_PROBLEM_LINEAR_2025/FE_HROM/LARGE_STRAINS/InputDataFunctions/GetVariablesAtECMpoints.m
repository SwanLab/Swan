function  MATPRO = GetVariablesAtECMpoints(MATPRO,fff,ngausT,ECMdata,setPointsElement,nTOTALdofsGAUSS,setIndices)



for i = 1:length(fff)
    NameProp = fff{i}  ;
    LENGTH_prop = size(MATPRO.(NameProp),1) ;
    if  LENGTH_prop == ngausT
        if isempty(ECMdata.setPoints)
            % MATPRO.(NameProp) = MATPRO.(NameProp)(ECMdata.setElements,:) ;
            % % CECM, ERROR FOUND  in 10-DEC-2021
            
            MATPRO.(NameProp) = MATPRO.(NameProp)(setPointsElement,:) ;  % CECM, change 10-Dec-2021
        else
            MATPRO.(NameProp) = MATPRO.(NameProp)(ECMdata.setPoints,:) ;  % DECM
        end
    elseif LENGTH_prop == nTOTALdofsGAUSS
        MATPRO.(NameProp) = MATPRO.(NameProp)(setIndices,:) ;
    else
        MATPRO.(NameProp) = MATPRO.(NameProp)(ECMdata.setElements,:) ;
    end
end