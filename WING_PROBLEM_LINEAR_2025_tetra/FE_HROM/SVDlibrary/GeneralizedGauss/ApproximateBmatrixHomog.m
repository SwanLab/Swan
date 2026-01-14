function [Bhred,BDh,IND_gausspointsCLOSE,DATAIN,coor_PARENTDOMAIN] =...
    ApproximateBmatrixHomog(COOR,CN,TypeElement,Ared,D,xNEW,nstrain,ngausELEM,VAR_SMOOTH_FE,DATAIN,DATAdistorsionELEM)

VAR_SMOOTH_FE = DefaultField(VAR_SMOOTH_FE,'ELEMENTS_xNEW',[]) ;
ELEMENTS_xNEW = VAR_SMOOTH_FE.ELEMENTS_xNEW ;
npoints = size(xNEW,1) ;
VAR_SMOOTH_FE = DefaultField(VAR_SMOOTH_FE,'B_MATRIX_FOR_COMPUTING_K_USING_INTERPOLATED_VALUE',0) ; % not used

% Loop over new points
nmodes = size(Ared,2) ;
nmodesD= size(D,2) ;
Bhred = zeros(nstrain*npoints,nmodes) ;
BDh = zeros(nstrain*npoints,nmodesD) ;

IND_gausspointsCLOSE = zeros(npoints,1) ;
errorBinterp = zeros(1,npoints) ; 
coor_PARENTDOMAIN =  zeros(size(xNEW)) ;

INDfin = 0 ;
for inew = 1:npoints
    xLOC = xNEW(inew,:);
    
    %     if isempty(ELEMENTS_xNEW)
    %
    %         INDnear = INDEXES_NEAR(inew) ;
    %         % Elements sharing INDnear.
    %         [ELEMnear, aaaa ]= find(CN == INDnear) ;
    %         ielem = 1;
    %         elemCONTAINER = [] ;
    %         % Searching for the element containing xLOC --> elemCONTAINER
    %         while ielem <= length(ELEMnear)
    %             elemLOC = ELEMnear(ielem) ;
    %             %         CNloc = CN(elemLOC,:);
    %             %         INDnodes =CNloc(IND_POLYG) ;
    %             %
    %             %         COORelem = COOR(INDnodes,:) ;
    %             %         %   COORelem = [COORelem; COORelem(1,:)] ;
    %             % %         plot(COORelem(:,1),COORelem(:,2))
    %             % %         plot(xLOC(1),xLOC(2),'*')
    %             %         [inELEM,onELEM] = inpolygon(xLOC(1),xLOC(2),COORelem(:,1),COORelem(:,2)) ;
    %             [inELEM,onELEM] = IsInsideGeneral(xLOC,COOR,CN,elemLOC,IND_POLYG)  ;
    %
    %             if inELEM == 1 || onELEM == 1
    %                 elemCONTAINER  = elemLOC ;
    %             end
    %             ielem = ielem + 1;
    %         end
    %
    %         if isempty(elemCONTAINER)
    %             error('Function for finding container elements did not work....')
    %         else
    %             IND_CONTAINING_ELEMENTS(inew ) = elemCONTAINER ;
    %         end
    %
    %     else
    elemCONTAINER = ELEMENTS_xNEW(inew) ;
    %  end
    
    % COMPUTATION OF THE LOCAL  B-MATRIX OF THIS POINT
    % -------------------------------------------------
    % -------------------------------------------------
    % Now we have the element defined by the nodes NODESloc = CN(elemCONTAINER,:)
    % with coordinates COOR(NODESloc,:), which contains xLOC
    indGAUSS = small2large(elemCONTAINER,ngausELEM) ;
    IND_gausspointsCLOSE(inew) = indGAUSS(1) ;
    NODESloc = CN(elemCONTAINER,:) ;
    
       METHOD_OLD =0;
    
    if METHOD_OLD == 1
        % Direct evaluation through quadratic functions
        [N,  BeTILDE] = ComputeBmatrixDirectly(COOR(NODESloc,:),xLOC,TypeElement) ;
    else
        DATAinp.ipoint = inew ;
        [N,  BeTILDE,coor_PARENTDOMAIN(inew,:)] = ShapeFunDer_inversemapping(COOR(NODESloc,:)',xLOC,TypeElement,DATAinp)  ;
        
    end
    
    
     
    
    
    
    ndim = size(COOR,2) ;
    % number of stain components
    
    Bpoint = QtransfB(BeTILDE,ndim) ;
    if size(Bpoint,1) ~= nstrain
        % Plane strain
        Bpoint = [Bpoint; zeros(1,size(Bpoint,2))] ;
    end
    
    % 
    
    % Contribution to the strain at this point
    
    DOFSloc = small2large(NODESloc,ndim) ;
    Bpoint_hyper = Bpoint*Ared(DOFSloc,:) ;
    BD_hyper = Bpoint*D(DOFSloc,:) ;
    INDINI = INDfin +1;
    INDfin = INDINI +nstrain-1;
    Bhred(INDINI:INDfin,:) = Bpoint_hyper ;
    BDh(INDINI:INDfin,:) = BD_hyper ;
    
    
    errorLOC =  Bpoint_hyper-VAR_SMOOTH_FE.Bhred_interp{inew} ;  
    errorLOC = norm(errorLOC,'fro') ; 
    errorBinterp(inew) =  errorLOC/norm(Bpoint_hyper)*100 ; 
end

% 
 DATAIN.MSGPRINT{end+1} = '---------------------------------------------' ; 
 disp( DATAIN.MSGPRINT{end} )
  DATAIN.MSGPRINT{end+1} = ['ELEMENT = ',num2str(ELEMENTS_xNEW(:)')] ; 
    DATAIN.MSGPRINT{end+1} = ['Distorsion angle (deg.) = ',num2str(DATAdistorsionELEM.distorsionANGLE(ELEMENTS_xNEW)')] ; 
  
  
   disp( DATAIN.MSGPRINT{end} )
    DATAIN.MSGPRINT{end+1} = ['Error Interp. Bhred (%) = ',num2str(errorBinterp(:)')] ; 
     disp( DATAIN.MSGPRINT{end} )

