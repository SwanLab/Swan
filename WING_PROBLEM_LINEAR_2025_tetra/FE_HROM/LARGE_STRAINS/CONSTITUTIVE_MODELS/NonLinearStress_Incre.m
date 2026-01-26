function VAR = NonLinearStress_Incre(DATA,OPERFE,VAR,FgradST)
%--------------------------------------------------------------------------
% NonLinearStress_Incre
%--------------------------------------------------------------------------
% PURPOSE
%   Extract the *nonlinear (incremental / inelastic)* part of the stress by
%   subtracting an *elastic reference response* computed with the **initial**
%   (virgin) elastic tangent `OPERFE.celastST_ini`. This is used when the
%   ECM/hyperreduction is applied **only to nonlinear stresses**.
%
% CONTEXT
%   • Set DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 1 upstream to activate.
%   • For small‑strain kinematics:   PK2_incre = PK2_total − C0 : E
%   • For large‑strain kinematics:   PK1_incre = PK1_total − C0 : (F − I)
%     (here, the “elastic” part is built with the initial tangent C0 and
%     kinematics consistent with your workflow).
%
% INPUTS
%   DATA      : structure with solver/config flags
%               - SMALL_STRAIN_KINEMATICS (1/0)
%               - CECM_ONLY_FOR_NONLINEAR_STRESSES (1/0)
%   OPERFE    : FE operators and elastic reference data
%               - celastST_ini : block matrix storing the initial C0
%                                (arranged by nstrain×nstrain blocks)
%               - Bst, IDENTITY_F (for large‑strain branch context)
%   VAR       : state struct with current fields already computed
%               - GLSTRAINS     : Green–Lagrange strains (small strain branch)
%               - PK2STRESS     : total PK2 stresses   (small strain branch)
%               - PK1STRESS     : total PK1 stresses   (large strain branch)
%               - DISP          : (used only in commented alternative)
%   FgradST   : deformation gradient at stress points (large‑strain branch)
%
% OUTPUTS
%   VAR.PK2STRESS_incre : (if small strain)  nonlinear (PK2) component
%   VAR.PK1STRESS_incre : (if large  strain) nonlinear (PK1) component
%
% IMPLEMENTATION NOTES
%   • celastST_ini is accessed “by blocks” using the stride pattern:
%       istrainCOL = istrain : nstrain : size(C0,1)
%     This treats the global (ngp*nstrain)-vectorized layout as stacked
%     blocks (one block per strain component).
%   • The commented alternatives show equivalent formulations using prebuilt
%     block‑diagonal products (ConvertBlockDiag) or Bst*DISP routes.
%   • For large strains there is a TODO about symmetrization (see note).
%
% REFERENCES (internal paths kept for traceability)
%   • 19_ExactLinearStiff.mlx  (small strain reference elastic part)
%   • MultiSnapStressFromDispNECMlarg.m (symmetrization remark)
%--------------------------------------------------------------------------

% Only act if we are in the regime where ECM targets nonlinear stresses.
if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1

    %==============================================================
    % SMALL‑STRAIN KINEMATICS BRANCH
    %==============================================================
    if DATA.SMALL_STRAIN_KINEMATICS == 1
        % Build the purely elastic PK2 response using the *initial* tangent C0,
        % then subtract it from the total PK2 to isolate the nonlinear part.
        %
        % Layout conventions:
        %   • celastST_ini is (nstrain*ngp) × (nstrain) “tall” per block column
        %     when accessed with the stride below.
        %   • VAR.GLSTRAINS is stacked by component: [E11; E22; E12; ...] per gp.

        % Pre-allocate elastic PK2 with same size as total PK2 vector
        STRESS_elastic = zeros(size(VAR.PK2STRESS));

        % nstrain is the block size (e.g., 3 for plane stress/strain, 6 in 3D)
        nstrain = size(OPERFE.celastST_ini,2);

        % Blocked multiplication: STRESS_elastic = C0 : GLSTRAINS
        for istrain = 1:nstrain
            % Indices of rows corresponding to output component 'istrain' across gps
            istrainCOL = istrain : nstrain : size(OPERFE.celastST_ini,1);

            for jstrain = 1:nstrain
                % Indices of columns corresponding to input component 'jstrain' across gps
                jstrainCOL = jstrain : nstrain : size(OPERFE.celastST_ini,1);

                % Accumulate C0(istrain,jstrain)*E_{jstrain} into σ_{istrain}
                STRESS_elastic(istrainCOL) = STRESS_elastic(istrainCOL) ...
                    + OPERFE.celastST_ini(istrainCOL,jstrain) .* VAR.GLSTRAINS(jstrainCOL);
            end
        end

        % Nonlinear (incremental) PK2 part to be hyperreduced/used by ECM
        VAR.PK2STRESS_incre = VAR.PK2STRESS - STRESS_elastic;

        % --- Alternative (kept for reference) ---------------------------------
        % Using a prebuilt block‑diagonal C0 and Bst*DISP:
        % celastST_iniD   = ConvertBlockDiag(OPERFE.celastST_ini);
        % OPERFE.celastINI_Bst = celastST_iniD * OPERFE.Bst;
        % VAR.PK2STRESS_incre  = VAR.PK2STRESS - OPERFE.celastINI_Bst * VAR.DISP;
        % ----------------------------------------------------------------------

    %==============================================================
    % LARGE‑STRAIN KINEMATICS BRANCH
    %==============================================================
    else
        % Build an “elastic” PK1 response from C0 applied to (F − I),
        % then subtract from the total PK1.
        %
        % NOTE (from your code): “Should we symmetrize this?” ECM works
        % with symmetrized gradients in some workflows. If your online stage
        % expects sym(∇u) or Green–Lagrange measures, consider adapting
        % (F − I) accordingly before production runs.

        % Pre-allocate elastic PK1 with same size as total PK1 vector
        STRESS_elastic = zeros(size(VAR.PK1STRESS));

        % nstrain equals the number of independent strain components per gp
        nstrain = size(OPERFE.celastST_ini,2);

        % Blocked multiplication: STRESS_elastic = C0 : (F − I)
        for istrain = 1:nstrain
            istrainCOL = istrain : nstrain : size(OPERFE.celastST_ini,1);

            for jstrain = 1:nstrain
                jstrainCOL = jstrain : nstrain : size(OPERFE.celastST_ini,1);

                % Accumulate C0(istrain,jstrain) * (F − I)_{jstrain}
                STRESS_elastic(istrainCOL) = STRESS_elastic(istrainCOL) ...
                    + OPERFE.celastST_ini(istrainCOL,jstrain) .* ...
                      (FgradST(jstrainCOL) - OPERFE.IDENTITY_F(jstrainCOL));
                % TODO (from original note): consider symmetrization if your ECM
                % derivation is based on symmetrized gradients.
            end
        end

        % Nonlinear (incremental) PK1 part for ECM / hyperreduction
        VAR.PK1STRESS_incre = VAR.PK1STRESS - STRESS_elastic;

        % --- Alternative (kept for reference) ---------------------------------
        % Using a converted block‑diag C0 acting on Bst*DISP:
        % celastST_iniD = ConvertBlockDiag(OPERFE.celastST_ini);
        % VAR.PK1STRESS_incre = VAR.PK1STRESS - celastST_iniD * (OPERFE.Bst * VAR.DISP);
        % ----------------------------------------------------------------------
    end
end

 

% Below is the original function, before introducing the in-line comments
% by ChatGPT-5 (23-August-2025)


% function VAR = NonLinearStress_Incre(DATA,OPERFE,VAR,FgradST)
% 
% if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1
%     
%     
%     if DATA.SMALL_STRAIN_KINEMATICS==1 ;
%         % THIS IS FOR SMALL STRAINS
%         % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
%         %   celastST_iniD = ConvertBlockDiag(OPERFE.celastST_ini) ; % Diagonal block matrix
%         % OPERFE.celastINI_Bst= celastST_iniD*OPERFE.Bst;
%         
%         %         OPTION = 0;
%         %
%         %         if OPTION == 1
%         STRESS_elastic = zeros(size(VAR.PK2STRESS)) ;
%         nstrain = size(OPERFE.celastST_ini,2) ;
%         for istrain = 1:nstrain
%             istrainCOL = istrain:nstrain:size(OPERFE.celastST_ini,1) ;
%             for jstrain = 1:nstrain
%                 jstrainCOL = jstrain:nstrain:size(OPERFE.celastST_ini,1) ;
%                 STRESS_elastic(istrainCOL) =   STRESS_elastic(istrainCOL) + OPERFE.celastST_ini(istrainCOL,jstrain).*VAR.GLSTRAINS(jstrainCOL) ;
%             end
%         end
%             VAR.PK2STRESS_incre = VAR.PK2STRESS - STRESS_elastic ;
% 
%         %         else
%         %
%         %             VAR.PK2STRESS_incre = VAR.PK2STRESS  -OPERFE.celastINI_Bst*VAR.DISP  ;
%         %         end
%         
%         
%     else
%         %   large strains
%         %
%         %         METHOD_loop = 1;
%         %
%         %         if METHOD_loop == 1
%         
%         STRESS_elastic = zeros(size(VAR.PK1STRESS)) ;
%         nstrain = size(OPERFE.celastST_ini,2) ;
%         for istrain = 1:nstrain
%             istrainCOL = istrain:nstrain:size(OPERFE.celastST_ini,1) ;
%             for jstrain = 1:nstrain
%                 jstrainCOL = jstrain:nstrain:size(OPERFE.celastST_ini,1) ;
%                 
%                 
%                 STRESS_elastic(istrainCOL) =   STRESS_elastic(istrainCOL) ...
%                     + OPERFE.celastST_ini(istrainCOL,jstrain).*(FgradST(jstrainCOL)-OPERFE.IDENTITY_F(jstrainCOL)) ;
%                 % 9-oct-2024...SHOULDN T WE SYMMETRIZED THIS? ECM WORKS
%                 % WITH SYMMETRIZED GRADIENTS, SEE
%                 % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/MultiSnapStressFromDispNECMlarg.m
%                 
%             end
%         end
%         
%         VAR.PK1STRESS_incre = VAR.PK1STRESS - STRESS_elastic ;
%         
%         %         else
%         %             % SEcond method
%         %             celastST_ini = ConvertBlockDiag(OPERFE.celastST_ini) ;
%         %
%         %             VAR.PK1STRESS_incre = VAR.PK1STRESS - celastST_ini*(OPERFE.Bst*VAR.DISP) ;
%         %
%         %         end
%         
%         
%         
%         
%         
%         
%         
%         
%     end
%     
%     
% end
