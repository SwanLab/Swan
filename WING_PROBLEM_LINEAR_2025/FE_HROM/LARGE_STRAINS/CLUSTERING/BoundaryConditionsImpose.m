function VAR = BoundaryConditionsImpose(VAR,DISP_CONDITIONS_cl,Fbody_cl,Ftrac_cl,iCL,istep)

dR = DISP_CONDITIONS_cl{iCL}.dR.U*DISP_CONDITIONS_cl{iCL}.dR.a(:,istep);
VAR.DISP(DISP_CONDITIONS_cl{iCL}.DOFr) = dR;

%if isstruct(Fbody_cl)
% 1.b) External forces
VAR.FEXT = Fbody_cl.U{iCL}*Fbody_cl.a(:,istep) + Ftrac_cl.U{iCL}*Ftrac_cl.a(:,istep) ;
    
%else

% 1.b) External forces
%VAR.FEXT = Fbody_cl{iCL}.U*Fbody_cl{iCL}.a(:,istep) + Ftrac_cl{iCL}.U*Ftrac_cl{iCL}.a(:,istep) ;
%end

 