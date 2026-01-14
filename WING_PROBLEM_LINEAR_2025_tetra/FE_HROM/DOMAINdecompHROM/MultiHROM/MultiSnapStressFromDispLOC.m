function [SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,VAR,ReactFORCE] = ...
    MultiSnapStressFromDispLOC(VAR,d,DATA,VARint_n,OTHER_output,OPERFE,MATPRO,iloc,...
    SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,iproj,INFO_RVE,Fbody_loc,Ftrac_loc)

if nargin == 0 
    load('tmp1.mat')
end


VAR.DISP  =d ;


if isempty(DATA.ListFieldInternalVariables)
    VARint_n = [] ;
else
    if iloc == 1
        for ivar = 1:length(DATA.ListFieldInternalVariables)
            FLOC = DATA.ListFieldInternalVariables{ivar};
            
            
            switch FLOC
                % Temporal amendment !! 
                case {'YieldStress','InternalVarStrain'}
                VARint_n.(FLOC) = OTHER_output.INICOND.(FLOC)(INFO_RVE.GaussINDEX_scalarSTR,:)  ;
             
                case {'PlasticStrains'}
               VARint_n.(FLOC) = OTHER_output.INICOND.(FLOC)(INFO_RVE.GaussINDEX_stress,:)  ;     
            end
                VAR.(FLOC) =  VARint_n.(FLOC)  ;
            
            
        end
    else
        for ivar = 1:length(DATA.ListFieldInternalVariables)
            FLOC = DATA.ListFieldInternalVariables{ivar};
            VARint_n.(FLOC) = VAR.(FLOC)  ;
        end
    end
    
    
end
DATA = DefaultField(DATA,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ; 
% STRESSES FROM HISTORY OF DISPLACEMENTS
[VAR,~,~,~] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
SNAPstressSTWOproj_LOC{iloc} = VAR.PK2STRESS ;
if  isempty(VAR.PK1STRESS)
    SNAPstressPonePROJ_LOC{iproj} = VAR.PK2STRESS ; ;
else
    SNAPstressPonePROJ_LOC{iproj} = VAR.PK1STRESS ;
end


% Internal forces (FOR ALL TIME STEPS)
% ---------------
VAR.FINT = InternalForces(OPERFE,VAR.PK1STRESS,VAR.PK2STRESS,DATA) ;

% ------------------------------------------------------------------------------
% REACTIVE FORCES (DIFFERENCE BETWEEN INTERNAL FORCES AND EXTERNAL FORCES)
% ------------------------------------------------------------------------------
% ONLY STATIC CASE SO FAR (8-Feb-2023)
% --------------------------------------------------
Fext = 0 ; 
if ~isempty(Fbody_loc.U)
    Fext = Fext + Fbody_loc.U(INFO_RVE.DOFS_globNUM,:)*Fbody_loc.a; 
end
if ~isempty(Ftrac_loc.U)
    Fext = Fext + Ftrac_loc.U(INFO_RVE.DOFS_globNUM,:)*Ftrac_loc.a; 
end 

ReactFORCE= VAR.FINT-Fext ; 

