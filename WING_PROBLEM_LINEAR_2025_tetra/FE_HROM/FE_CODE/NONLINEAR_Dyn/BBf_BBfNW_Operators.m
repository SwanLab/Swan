function [BBf,BBfNW,massMf,massMs] = BBf_BBfNW_Operators(OPERfe,BBnw,BB,massM) 

  
 
 
if isempty(OPERfe.Gbound)
    OPERfe.Gbound = 0 ;
end
if ~isempty(OPERfe.DOFm)
    if ~isempty(BB)
    BBf = [BB(:,OPERfe.DOFl)   BB(:,OPERfe.DOFm)+BB(:,OPERfe.DOFs)*OPERfe.Gbound] ;
    else
       BBf = [] ;  
    end
    BBfNW = [BBnw(:,OPERfe.DOFl)   BBnw(:,OPERfe.DOFm)+BBnw(:,OPERfe.DOFs)*OPERfe.Gbound] ;
%     if ~isempty(ASSEMBLY_INFO)
%         error('Assembly method only valid when DATA.AssemblyMethodK_BCB = 1')
%     end
else
    if ~isempty(BB)
    BBf = BB(:,OPERfe.DOFl)  ;
    
    else
        BBf = [] ; 
     end
    BBfNW = BBnw(:,OPERfe.DOFl)  ;
   
  %  if ~isempty(ASSEMBLY_INFO)
   %     ASSEMBLY_INFO.DOFl = OPERfe.DOFl ; 
   % end
end
if ~isempty(massM)
    [massMf massMs] = MassBlockMatrices(OPERfe.DOFl,OPERfe.DOFm,OPERfe.DOFs,OPERfe.Gbound,massM) ;
    % Acceleration of BCs
    %   nDOFTOT = size(COOR,1)*size(COOR,2);
    %   [d0,v0,gBOUNDdd,gBOUNDd]=AccelDirichBC(DynDATA,nDOFTOT,DOFm,DOFs,DOFl,Gbound,gBOUND,VECTOR_TIME);
    %else
    %  d0 = [] ; v0=[] ; gBOUNDdd = [] ;
else
    massMf =[] ;
    massMs =[] ;
end