function setelem = SetElementsBoundary(CONNECTb,CNb)


METHOD_CONNECTIVITIES_VECTORISED = 1; % Default method.
%
if METHOD_CONNECTIVITIES_VECTORISED == 1
    
    % Problem detected 7th-January-2019
    % This function does not work if the connectivity matrix has repeated
    % elements, hence the below conditional
    % Nevertheless, in MeshGenerationRepeat3D_faces.m, we have included
    % a snippet of code that detects repeated elements. Thus, the
    % conditional below may be redundant.
    
    % ------------------------------
    [dummy setelem]= ElemBnd(CONNECTb,unique(CNb(:)));
    if length(setelem) ~=size(CNb,1)
        % AMendment 7-Jan-2019
        CnbNEW = CONNECTb(setelem,:) ;
        oldNODES = unique(CNb(:));
        newNODES = unique(CnbNEW(:));
        warning('ERROR in preparing GID data')
        disp('There is something wrong in the construction of CNb')
        disp('Perhaps non-interface forces are ill-defined')
        disp('Revise the definition of faces in the GID-preprocess file')
        disp('If there is nothing wrong there, then turn METHOD_CONNECTIVITIES_VECTORISED =0')
        error(' ')
    end
    
    % August -2019
    % Now we have to order
    CNbNEW = CONNECTb(setelem,:) ;
    % so that it maches CNb. In principle, we can do this by
    Locb1 = {};   % Change 4-Feb-2020
    Lia = {} ;
    for  iiii = 1:size(CNb,2)
        [Lia{iiii},Locb1{iiii}] = ismember(CNb(:,iiii),CNbNEW(:,iiii)) ;
    end
    
    salir =0 ;
    iiii=1 ;
    while salir == 0 && iiii < size(CNb,2)
        
        setelemTRY = setelem(Locb1{iiii}) ;
        if sum(sum(abs(CONNECTb(setelemTRY,:)-CNb))) == 0  && iiii <=size(CNb,2)
            setelem = setelemTRY ;
            salir = 1;
        else
            iiii = iiii+1;
        end
        
        
        
        
    end
    
    if iiii > size(CNb,2)
        error('')
    end
    
else
    warning('This routine requires a vectorized version !!!!!')
    % Non-vectorized version. Very slow for large number of elements
    % It was employed in favor of the vectorized version
    % from July 2016 till January 2019 (because we didn't realize that there
    % was an issue with repeated elements in CONNECTb)
    
    % It was again enabled on 23-June-2019 because we realized that the
    % vectorized version altered the results.
    
    % 23-Aug-2019. The confusion lies in that  ElemBndSort and ElemBnd
    % perform totally different operations. We are interested in the
    % operation performed by ElemenBndSort, not by ElemBnd. Indeed, our
    % goal is to determine which are the entries of CONNECTb corresponding
    % to CNb ---this is what is done in ElemeBndSort. By contrast, ElemBnd
    % provides the entries of CONNECTb corresponding to a given set of
    % indices.
    
    % Therefore, all we can do is to VECTORIZED ElemBnd ---> This is done by adding
    
    %  CNbNEW = CONNECTb(setelem,:) ;
    % so that it maches CNb. In principle, we can do this by
    % [Lia,Locb] = ismember(CNb(:,1),CNbNEW(:,1)) ;
    % setelem = setelem(Locb) ;
    %
    
    % in method 1.
    
    [ setelem] = ElemBndSort(CONNECTb,CNb) ;    % CNb is a subset of CONNECTb.
    
    
end