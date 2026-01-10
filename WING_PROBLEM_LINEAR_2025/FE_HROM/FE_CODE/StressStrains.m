function [strainGLO stressGLO posgp DATAOUT ]= StressStrains(COOR,CN,TypeElement,celasglo...
    ,d,typePROBLEM,celasgloINV,Bst,Cglo,DATA,wST,Nst,React)
%%%% COmputation of strain and stresses   at each gauss point
%dbstop('5')
if nargin == 0
    load('tmp1.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
nstrain = DATA.nstrain ;
DATA = DefaultField(DATA,'posgp_given',[]) ; 
% 4-Nov-2022 --- GAuss points given provided by the user (see for instance
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/MultiscaleHROM/linRECT_CROSS_SECTION/DATAFE_slice_linear3D.m
%
if isempty(DATA.posgp_given)
TypeIntegrand = 'K' ;
else
   TypeIntegrand = {DATA.posgp_given,DATA.weights_given} ;  
end
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;


ngaus = length(weig) ;

DATAOUT.stressVONMISES = [] ;

DATAOUT.stressAVG = [] ;
DATAOUT.weig = weig ;

DATA  =DefaultField(DATA,'PRINT_AVERAGE_STRESSES_ON_ELEMENTS',1) ;

%dbstop('22')
if  DATA.BCBformulation==1
    
    [strainGLOv,stressGLOv,DATAOUT,strainGLO,stressGLO,ngausLOC,wST] = ...
        StressStrains_vect(d,Bst,Cglo,DATA,wST,nelem,nstrain,ngaus,ndim) ;
else
    
    % Non-vector format
    % --------------------------------
    switch typePROBLEM
        case 'pstrain'
            TP = 1;
            %    if nstrain == 3
            nstrain=4 ;
            %   end
        case 'pstress'
            TP = 2 ;
        case '3D'
            TP = 3 ;
    end
    
    stressGLO = zeros(ngaus*nstrain,nelem);
    strainGLO = zeros(ngaus*nstrain,nelem);
    for e = 1:nelem
        % Elasticity matrix of element "e"
        celas = celasglo(:,:,e) ;
        celas3Dinv = celasgloINV(:,:,e) ;
        celas3D = inv(celas3Dinv) ;
        % Coordinates of the nodes of element "e"
        CNlocNOD = CN(e,:) ;
        Xe = COOR(CNlocNOD,:)' ;
        % Displacement at   nodes of element e
        CNloc = Nod2DOF(CNlocNOD,ndim) ;
        dE = d(CNloc) ;
        
        
        
        
        for  g = 1:ngaus
            % Matrix of derivatives for Gauss point "g"
            BeXi = dershapef(:,:,g) ;
            % Jacobian Matrix
            Je = Xe*BeXi' ;
            % Matrix of derivatives with respect to physical coordinates
            BeTILDE = inv(Je)'*BeXi ;
            % Matrix of symmetric gradient
            Be = QtransfB(BeTILDE,ndim) ;
            %
            strain = Be*dE;
            %  stress = celas*strain ;
            % Additional component
            if TP == 1
                strain3D = [strain(1) strain(2) 0 0 0 strain(3)]'; % Plane strain
                stress3D = celas3D*strain3D ;
                stressGID = stress3D([1 2 6 3]) ;  % For post-processing with GID
                strainGID = strain3D([1 2 6 3]) ;
            elseif TP == 2
                stress = celas*strain ;
                stress3D = [stress(1) stress(2) 0 0 0 stress(3)]'; % Plane strain
                strain3D = celas3Dinv*stress3D ;
                
                strainGID = strain3D([1 2 6 3]) ;
                stressGID = stress3D([1 2 6 3]) ; % For post-processing with GID
            elseif TP ==3
                stress = celas*strain ;
                
                stressGID = stress([1 2 3 6 4 5]) ; % For post-processing with GID
                strainGID = strain([1 2 3 6 4 5]) ; % For post-processing with GID
            end
            
            
            indINI = (g-1)*nstrain+1 ;  indFIN = nstrain*g;
            stressGLO(indINI:indFIN,e) =  stressGID ;
            strainGLO(indINI:indFIN,e) =  strainGID ;
        end
    end
end


% Average stresses and strains
%dbstop('120')
if DATA.CALCULATE_averageSTRESS==1
    % 1st order homogenization
    DATAOUT = AverageStressStrain_1st(DATA,stressGLOv,wST,nstrain,ngausLOC,nelem,DATAOUT,strainGLOv) ;
elseif DATA.CALCULATE_averageSTRESS==2
    DATAOUT = AverageStressStrain_2nd(DATA,stressGLOv,wST,nstrain,ngausLOC,nelem,DATAOUT,strainGLOv,...
        Nst,TypeElement,nnodeE,ndim,COOR,CN,ngaus) ;
    
elseif DATA.CALCULATE_averageSTRESS==3
    
    if ~isempty(React)
        
        DATAOUT = AverageStressGEN_react(DATA,COOR,React,DATAOUT) ;
        
    else
        error('Option not compatible')
        
    end
    
end