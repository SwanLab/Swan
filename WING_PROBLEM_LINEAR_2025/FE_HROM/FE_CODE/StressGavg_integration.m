function [stressAVG Qx_avg Qy_avg VOLUME]= StressGavg_integration(DATA,stressGLOv,wST,nstrain,ngaus,nelem,DATAOUT,strainGLOv,...
    Nst,TypeElement,nnodeE,ndim,COOR,CN,ngausLOC)

%dbstop('5')
if nargin == 0
    load('tmp.mat')
end

Qx_avg = [] ; Qy_avg = [] ;
stressGLOv_W = bsxfun(@times,stressGLOv,wST) ;
stressGLOv_W =reshape(stressGLOv_W,nstrain,[]) ;

strainGLOv_W =bsxfun(@times,strainGLOv,wST) ;
strainGLOv_W =reshape(strainGLOv_W,nstrain,[]) ;
strainGLOv =reshape(strainGLOv,nstrain,[]) ;


stressAVG = sum(stressGLOv_W,2) ;
% Strains

% Computing volume

if DATA.NOVOIDS == 1 ;
    ROWS = 1:nstrain:nstrain*ngausLOC*nelem ;
    wSTs = wST(ROWS);
    VOLUME = sum(wSTs) ;
else
    error('Option not implemented')
end

stressAVG_comp = stressAVG/VOLUME ;



%%%5 Moments
% We need the z-coordinates of all gauss points

if  DATA.RECALCULATE_STIFFNESS == 1
    
    if isempty(Nst)
        % Calculating Nst matrix
        [ Nelem  ] = ComputeNelemALL(TypeElement,nnodeE,ndim,nelem) ;
        % -------------------------------------------------------------------------
        disp('Done')
        % Assembly of Bst
        nnode = size(COOR,1) ;
        disp('Assembly of Nst (stacked shape functions matrix)...')
        Nst = AssemblyNGlobal(Nelem,nelem,nnodeE,ndim,ngaus,CN,nnode) ;
        disp('Done')
        
    end
    if DATA.STORE_STIFFNESS ==1
        disp('Storing Nst Matrix...')
        save(DATA.nameWORKSPACE,'Nst','-append');
        disp('Done')
    end
    
else
    load(DATA.nameWORKSPACE,'Nst')
    
end






%%%%

%dbstop('44')
X = COOR' ; X = X(:);
Xgauss = Nst*X ;
% z- coordinates
Z = Xgauss(3:3:end);

%%% Computation of Mx, My, Mxy

% STRESSES = {stressGLOv(1,:)',stressGLOv(2,:),stressGLOv(6,:)};
% Mx = int(sigma_x,z) divided by VOL, and multiplied by thickness
zmax = max(COOR(:,3)) ; zmax =zmax(1) ;
zmin = min(COOR(:,3)) ; zmin =zmin(1) ;
dZ = zmax-zmin ; A = VOLUME/dZ ;
Mx = sum(stressGLOv_W(1,:).*Z')/A ;
My = sum(stressGLOv_W(2,:).*Z')/A ;
Mxy = sum(stressGLOv_W(6,:).*Z')/A ;

%%%%

stressAVG = [stressAVG_comp(1)*dZ
    stressAVG_comp(2)*dZ
    stressAVG_comp(6)*dZ
    Mx
    My
    Mxy
    stressAVG_comp(4)*dZ
    stressAVG_comp(5)*dZ
    ];

%dbstop('100')
if DATA.ORDER_SHEAR_THEORY == 3  % ; DATA.SPECIAL_BC_SHEAR_ON == 3 | DATA.SPECIAL_BC_SHEAR_ON == 13 ;
    
    CoeffZH = 3/2*(1 - (4/dZ^2)*(Z.^2)') ;
    Qx = sum(stressGLOv_W(5,:).*CoeffZH)/A ;
    Qy = sum(stressGLOv_W(4,:).*CoeffZH)/A ;
    stressAVG(7) = Qy ;
    stressAVG(8) = Qx ;
    
    Qx_avg = stressAVG_comp(5)*dZ ;
    Qy_avg = stressAVG_comp(4)*dZ ;
    
%     PRUEBAS =1 ;
%    % dbstop('113')
%     if PRUEBAS==1
%         enerYZ = stressGLOv_W(4,:).*strainGLOv(4,:) ;
%         enerYZ = sum(enerYZ) ;
%         gamma_yz = DATA.strainINP(7) ;
%         supRVE = VOLUME/dZ ;
%         QyENER = enerYZ/gamma_yz/supRVE ;        
%       %s  stressAVG(7) = Qy ;
%     end
    
%else
 %   stressAVG = ShearCoefficient(DATA,stressAVG,stressGLOv_W,strainGLOv,VOLUME,dZ) ;
end
