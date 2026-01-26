function [UnrestrainedDOFS,a] = UnrestrainedDOFSlocal(a,iface,DOFB,ndim)

UnrestrainedDOFS =[] ; 
if iscell(a{iface})    
    NotConstrainedDOFS = cellfun(@isempty,a{iface}) ;  % Not constrained DOFs  
    if all(NotConstrainedDOFS)
        a{iface} = [] ;
    elseif all(NotConstrainedDOFS == 0)
        a{iface} = cell2mat(a{iface}) ;
    elseif  any(NotConstrainedDOFS) 
        IND_nc =  find(NotConstrainedDOFS == 1) ;  % Not contrained indices        
        if any(find(IND_nc > ndim))
            error('Rotations cannot be unrestricted')
        end        
        UnrestrainedDOFS = [] ; 
        for iii = 1:length(IND_nc)
            UnrestrainedDOFS = [UnrestrainedDOFS  ;  IND_nc(iii):ndim:length(DOFB)] ; 
            a{iface}{IND_nc(iii)} = 0 ; 
        end        
        a{iface} = cell2mat(a{iface}) ; 
    end
end