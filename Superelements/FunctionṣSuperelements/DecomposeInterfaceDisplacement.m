function [uI] = DecomposeInterfaceDisplacement(GlobalU,interfaces,subD,initialDof)
    uI = cell(1,interfaces); 
    for i=2:interfaces
        DofsInterface = (1:1:size(subD{i,1},2)) + initialDof;
        subuI = GlobalU(DofsInterface);
        subuI = [subuI(1:2:end-1), subuI(2:2:end)];
        uI(i) = {subuI};
        initialDof = max(DofsInterface);
    end
end