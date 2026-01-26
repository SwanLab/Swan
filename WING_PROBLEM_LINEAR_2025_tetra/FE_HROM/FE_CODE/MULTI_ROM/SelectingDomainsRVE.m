function [SELECTED_DOMAINS] =  ...
    SelectingDomainsRVE(DATAIN,ISJOINT,ientity,MAXstressVONMISES,NO_STRESSES,...
    SELECTED_DOMAINS,ELEMS)

nmax = length(ELEMS);

if ~isempty(DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER) & DATAIN.DO_NOT_COMPUTE_STRESSES ==0  & ISJOINT ==0
    % error('Modify this part of the code  !!!!!')
    NUMBERprint = min(DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER,nmax) ;
    switch DATAIN.DOMAINS_POSTPROCESS_SELECT.VARIABLE
        case 'VONMISES'
            [~,INDEXES] =  sort(MAXstressVONMISES,'descend') ;
            NUMBERprint = min(NUMBERprint,length(INDEXES)) ;
            SELECTED_DOMAINS_loc = unique(INDEXES(1:NUMBERprint)) ;
            SELECTED_DOMAINS_loc = [1,SELECTED_DOMAINS_loc,length(MAXstressVONMISES)] ; % First and last domains are included
            SELECTED_DOMAINS_loc = unique(SELECTED_DOMAINS_loc) ;
            SELECTED_DOMAINS{ientity} = SELECTED_DOMAINS_loc ;
            
            %             if  ~isempty( stressVONMISES_glo{ientity} )
            %             stressVONMISES_glo{ientity}  =  stressVONMISES_glo{ientity}(:,SELECTED_DOMAINS_loc) ;
            %             stressGLO_glo{ientity} =  stressGLO_glo{ientity}(:,SELECTED_DOMAINS_loc) ;
            %             end
            %      DISP3D = DISP3D(:,SELECTED_DOMAINS) ;
            %    DATAIN.DOMAINS_POSTPROCESS = SELECTED_DOMAINS ;
        otherwise
            error('Option not implemented')
    end
else
    if ~isempty(DATAIN.DOMAINS_POSTPROCESS) && NO_STRESSES == 0 && DATAIN.DO_NOT_COMPUTE_STRESSES ==0
        
        [DDD iA iB]= intersect(DATAIN.DOMAINS_POSTPROCESS,ELEMS) ;
        
        % stressVONMISES_glo{ientity}  =  stressVONMISES_glo{ientity}(:,iB) ;
        % stressGLO_glo{ientity}  =  stressGLO_glo{ientity}(:,iB) ;
        SELECTED_DOMAINS{ientity} = iB ;
    elseif  isempty(DATAIN.DOMAINS_POSTPROCESS) && NO_STRESSES == 1
        % stressVONMISES =  [] ;
        % stressGLO =  [] ;
        SELECTED_DOMAINS = [1:length(ELEMS)];
    end
    
    
end

end
