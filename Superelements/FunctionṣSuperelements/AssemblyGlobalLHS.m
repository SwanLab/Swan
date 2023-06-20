function [GlobalLHS] = AssemblyGlobalLHS(AGinputs)
    type = AGinputs.type;
    K = AGinputs.K;
    C = AGinputs.C;
    
    if type == "Two"
        GlobalLHS = [  K                   C           ;
                       C'   sparse(size(C,2),size(C,2))];
    elseif type == "Three"
        D = AGinputs.D;
        GlobalLHS = [                 K                               C                    sparse(size(C,1),size(D,2));
                                      C'                 sparse(size(C,2),size(C,2))                     D            ;
                      sparse(size(D,2),size(C,1))                     D'                   sparse(size(D,2),size(D,2))];
    end
end