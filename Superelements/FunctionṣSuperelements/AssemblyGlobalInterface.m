function [GlobalD] = AssemblyGlobalInterface(ADinputs)
    subD = ADinputs.subD;

    interfaces = size(subD,1);
    GlobalD = sparse(size(subD{2,2},1),size(subD{2,2},2));
    for i=2:interfaces
        GlobalD = [               GlobalD                              sparse(size(GlobalD,1),size(subD{i,1},2));
                   sparse(size(subD{i,1},1),size(GlobalD,2))                             -subD{i,1}             ;
                   sparse(size(subD{i,2},1),size(GlobalD,2))                              subD{i,2}             ];
    end
end