function     [stressSTNP1,Unp1,VELnp1,ACELnp1,EPnp1,alphaNP1,sigmayNP1,iter, CONVERGENCE,TIMEelap] =...
    NewtonRaphsonDyn_FE2D(tolNEWTONRAPSHON,EPn,sigmayN,alphaN,...
    BBfNW,istep ,Un,max_iter,F_f,DOFf,B_s_g0,PROPMAT,...
    massMf,massMs,gDD,DynDATA,VELn,ACELn,increT,BBfT,wST,ASSEMBLY_INFO)

% if istep == 49
% dbstop('9')
%  end
if nargin == 0
    load('tmp.mat')
end

% NEWTON-RAPHSON
norm_res = 1;
tol_rel  = tolNEWTONRAPSHON;
iter = 1 ;
SALIR = 0 ;


DynDATA = DefaultField(DynDATA,'betaNM',1/4) ;
DynDATA = DefaultField(DynDATA,'gammaNM',1/2) ;


betaNM  = DynDATA.betaNM ;
gammaNM = DynDATA.gammaNM;

% Initial displacements (for iteration purposes)
Unp1 = Un ;
if  DynDATA.DYNAMIC == 0
    U_k = Un(DOFf) ;
else
    uTIL  = Un(DOFf) + VELn*increT +  0.5*(1-2*betaNM)*ACELn*increT^2 ;
    U_k    =  uTIL ; %
end
deltaU = zeros(size(U_k)) ;

% if istep ==45
%     dbstop('34')
%     gammaNM
% end
Celas = [] ; % Elasticity matrix

TIMEelap.UpdateStress = 0 ;
TIMEelap.Residual = 0 ;
TIMEelap.AssesmblyK = 0 ;
TIMEelap.Solver = 0 ;
TIMEelap.Total = 0 ;

TimeTotal = tic;
while ((SALIR == 0) && (iter <= max_iter))
    
    
    
    % dbstop('39')
    eNP1 = BBfNW*U_k + B_s_g0  ; % Computing "strain-like" global vector field
    % Computing vector of stresses and tangent moduli at each gauss point
    % ------------------------------------------------------------------------
    
    % dbstop('39')
    CALC_CTANG = 1 ;
    
    TIMEloc = tic;
    if isempty(BBfT)
        WEIGHTS = wST ;
    else
        WEIGHTS = [] ;
    end
    
    
    [Ctang,   stressSTNP1, EPnp1,  alphaNP1,  sigmayNP1,CtangM] = ...
        UpdateStresses_J2_FE2D(eNP1,PROPMAT,EPn,sigmayN,alphaN,CALC_CTANG,iter,WEIGHTS) ;
    TIMEloc = toc(TIMEloc) ;
    TIMEelap.UpdateStress = TIMEelap.UpdateStress + TIMEloc ;
    
    
    
    %     if STORE_CBMAT_SNAPSHOTS == 2
    %         ifin = length(U_k)*iter; iini = length(U_k)*(iter-1)+1 ;
    %         CBMAT(:,iini:ifin) = CtangNP1_B ;
    %     end
    
    % dbstop('52')
    
    TIMEloc = tic;
    % ----------------
    if ~isempty(BBfT)
        FINT = BBfT*stressSTNP1 ;
    else
        FINT = BBfNW'*(stressSTNP1.*wST) ;
    end
    TIMEloc = toc(TIMEloc) ;
    
    
    
    %%%%
    % External forces
    if DynDATA.DYNAMIC == 0
        FEXT = F_f ;
        ACELnp1 = ACELn ;
    else
        % Dynamic case
        % ------------
        ACELnp1 =  (U_k-uTIL)/betaNM/increT^2;
        FEXT = F_f -massMs*gDD- massMf*ACELnp1 ;
    end
    %%%%
    %     if ~isempty(DATA_INTERPOLATORY)
    %         if DATA_INTERPOLATORY.BLOCKDECOMP_EXPMAT == 1
    %             FEXT = DATA_INTERPOLATORY.Dfext*FEXT ;
    %         end
    %     end
    % Residual
    ResF= FINT -FEXT;%  F_f;      % Computing residual vector
    %disp(['ResF (borrar luego) =',num2str(norm(ResF))]); 
    criter_f =  CheckConvergence(FINT,FEXT,ResF,num2str(iter));
    
    TIMEelap.Residual = TIMEelap.Residual + TIMEloc ;
    
    
    
    % -------------------------------------
    if criter_f<tol_rel  | norm_res <tol_rel
        SALIR = 1 ;
        U_kp1 = U_k ;
        %         % dbstop('78')
        %         [ResF,FINT_Z] = CalF_RESmethod (ISRESID,BBfNWTnp,stressSTNP1,F_fNP,massMsNP,gDD,massMfNP,ACELnp1,...
        %    DATA_FE1         nmodesU,CUBATURE_FINT_FEXT,FINT_Z,BBfNWT,FEXT,FINT) ;
    else
        
        
        TIMEloc = tic;
        
        if isempty(ASSEMBLY_INFO)
            
            IMPLE = 0;
            
            if IMPLE  == 1
                % 25-June-2019 ---- Attempt to improve the efficiency ....
                Kstiff = BBfNW'*Ctang*BBfNW ;
                %             NPARTS = 100 ;
                %             INDgauss = unique(ceil(linspace(1,size(BBfNW,1),NPARTS))) ;
                %             Kstiff = sparse(size(BBfNW,2),size(BBfNW,2)) ;
                %
                %             for iparts = 1:NPARTS-1
                %                 IND = INDgauss(iparts):INDgauss(iparts+1) ;
                %                 Kstiff = Kstiff + BBfNW(IND,:)'*Ctang(IND,IND)*BBfNW(IND,:) ;
                %             end
            else
                
                % dbstop('92')
                if ~isempty(BBfT)
                    CtangNP1_B = Ctang*BBfNW ;
                    Kstiff = BBfT*CtangNP1_B ; % Jacobian matrix
                else
                    Kstiff = Ctang*BBfNW ;
                    Kstiff = BBfNW'*Kstiff ;     % Jacobian matrix  (Weights included in CtangNP1_B)
                end
            end
        else
            Kstiff = AssemblyKMethodStandard_nonl(ASSEMBLY_INFO,CtangM) ;
        end
        
        
        
        if DynDATA.DYNAMIC ==1
            Kstiff = Kstiff + massMf/increT^2/betaNM  ;
        end
        
        TIMEloc = toc(TIMEloc) ;
        TIMEelap.AssesmblyK = TIMEelap.AssesmblyK  + TIMEloc;
        
        TIMEloc =tic ;
        dirSEARCH = - Kstiff\ResF ;
        TIMEloc = toc(TIMEloc) ;
        TIMEelap.Solver = TIMEelap.Solver + TIMEloc ;
        
        % Energ = -dirSEARCH'*ResF;
        
        % if ~isempty(DATA_INTERPOLATORY) &&  ~isempty(DATA_INTERPOLATORY.ALPHA_PERTURB)
        
        %             if Energ <0
        %                 factorr = norm(Kstiff);
        %                 alpha_per = DATA_INTERPOLATORY.ALPHA_PERTURB ;
        %                 Kstiff = (alpha_per)*factorr*eye(size(Kstiff)) + (1-alpha_per)*Kstiff ;
        %                 dirSEARCH = - Kstiff\ResF ;
        %                 Energ = -dirSEARCH'*ResF;
        %
        %             end
        %
        
        
        %     end
        
        
        U_kp1 = U_k + dirSEARCH;  % Displacement update
        %dbstop('136')
        %   deltaU =  U_kp1-U_k ;
        iter = iter+1 ;
        U_k = U_kp1 ;
        
        
    end
end



% TESTSOLUTION =1 ;
% if TESTSOLUTION ==1
%
%
%     CtangNP1_B = Ctang*BBfNW ;
%              dbstop('143')
%             Kstiff = BBfT*CtangNP1_B ;
%
%     CBs = BBfT*Ctang*B_s_g0
%     dU = -Kstiff\CBs ;
% end





%%% Updating displacements and velocities
Unp1(DOFf) = U_kp1 ;

if DynDATA.DYNAMIC == 1
    VEL_TILDE = VELn  + (1-gammaNM)*increT*ACELn;
    VELnp1 =  VEL_TILDE + gammaNM*increT*ACELnp1;
else
    
    VELnp1 = VELn ;
end

TimeTotal = toc(TimeTotal);

TIMEelap.Total = TimeTotal ;


%  if istep == 965
%         dbstop('158')
%         Kstiff;
%
%
% end


% %%%
%   DETERMINE_EIGENFREQ = 1 ;
% % %
% if DETERMINE_EIGENFREQ == 1
%   %  if istep == 8000
%         dbstop('173')
%         disp('Eigenvalues')
%         K =   (BBfT*CtangNP1_B )  ;
%        % M =  (massMf) ;
%         [DIAGONALMval] = eigs(massMf,K);
%         DIAGONALMval = 1./DIAGONALMval ;
%         EIGENFREQ = sqrt((DIAGONALMval)) ;
%         EIGENPERIOD = 2*pi./EIGENFREQ ;
%         EIGENPERIOD = sort(EIGENPERIOD,'descend') ;
%         EIGENPERIOD
%
%
%    % end
% end

%dbstop('186')
% if STORE_CBMAT_SNAPSHOTS == 1
%     CBMAT = CtangNP1_B ;
% elseif STORE_CBMAT_SNAPSHOTS ==3
%     CBMAT = CtangNP1_B*deltaU ;
% end

CONVERGENCE =1 ;

if iter> max_iter
    % dbstop('202')
    warning('convergence error')
    CONVERGENCE = 0 ;
    
end