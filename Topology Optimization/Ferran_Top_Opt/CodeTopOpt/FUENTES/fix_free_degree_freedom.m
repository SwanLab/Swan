function [fix_df,free_df]=fix_free_degree_freedom(type,nndof,nunkn,fixnodes,problembsc)

%Identifica y modifica grados de libertad y restringidos
    switch type
        case {'TRIANGLE','QUAD'}
            if (size(fixnodes,1)>0)
                fix_df = (fixnodes(:,1)-1)*nunkn + fixnodes(:,2);  %Finds the equation number
                free_df = setdiff ( 1:nndof, fix_df );
            else
                free_df = (1:nndof);
                fix_df = [];
            end   
        case 'LINEAR_TRIANGLE_MIX'
            if (size(fixnodes,1)>0)
                fix_df = (fixnodes(:,1)-1)*3 + fixnodes(:,2);  %Finds the equation number
                free_df = setdiff ( 1:nndof, fix_df );
            else
                free_df = (1:nndof);
                fix_df = [];
            end
        case 'LINEAR_TRIANGLE_MIX_COUPLED'
            fix_df = (fixnodes(:,1)-1)*3 + fixnodes(:,2);  %Finds the equation number
            free_df = setdiff ( 1:nndof, fix_df );
            
        case 'HEXAHEDRA'
            if (size(fixnodes,1)>0)
                fix_df = (fixnodes(:,1)-1)*3 + fixnodes(:,2);  %Finds the equation number
                free_df = setdiff ( 1:nndof, fix_df );
            else
                free_df = (1:nndof);
                fix_df = [];
            end
        otherwise
            error('No existe es tipo de elemento o no ha sido implementado')
    end

end