classdef ElasticProblemFluc < ElasticProblem

    properties (Access = public)
        sizePer
        variables2print
    end

    methods (Access = public)
        
        function Ch = computeLagrangeMultSum(obj)
            nstre = obj.material.nstre;
            basis = diag(ones(nstre,1));
            Ch = zeros(nstre,nstre);
            for istre=1:nstre
                obj.vstrain = basis(istre,:);
                obj.solve();

                perDOFslave = obj.boundaryConditions.periodic_constrained;
                obj.sizePer = size(perDOFslave, 1);
                nEqperType  = obj.sizePer/4;
                L           = obj.variables.LangMult;
                sigmaX      = 0;
                sigmaY      = 0;
                tauXY       = 0;
                d1          = obj.sizePer;
                LPer        = L(1:d1);
                LDir        = L(d1+1:end);
                for i = 1:nEqperType
                    sigmaX = sigmaX + LPer(i);
                end
                for i = nEqperType+1:2*nEqperType
                    tauXY  = tauXY + LPer(i);
                end
                for i = 2*nEqperType+1:3*nEqperType
                    tauXY  = tauXY + LPer(i);
                end
                for i = 3*nEqperType+1:4*nEqperType
                    sigmaY = sigmaY + LPer(i);
                end
                sigmaX       = sigmaX + LDir(1) + LDir(5);
                sigmaY       = sigmaY + LDir(2) + LDir(4);
                tauXY        = (tauXY + LDir(6) + LDir(3) - LDir(7) - LDir(8))/2;
                Ch(istre, :) = [sigmaX; sigmaY; tauXY];
                obj.variables.Chomog = -Ch;
            end
        end

        function Ch = computeChomog(obj)
             nelem = size(obj.material.C,3);
            npnod = size(obj.displacementFun.fValues,1);
            ndofs = npnod*obj.displacementFun.ndimf;
            nstre = obj.material.nstre;
            ngaus = obj.quadrature.ngaus;
            basis = diag(ones(nstre,1));
            tStrn  = zeros(nstre,ngaus,nstre,nelem);
            tStrss = zeros(nstre,ngaus,nstre,nelem);
            tDisp  = zeros(nstre,ndofs);
            Ch = zeros(nstre,nstre);
            for istre=1:nstre
                obj.vstrain = basis(istre,:);
                obj.solve();
                vars = obj.computeStressStrainAndCh();
                Ch(:,istre)         = vars.stress_homog;
                tStrn(istre,:,:,:)  = vars.strain;
                tStrss(istre,:,:,:) = vars.stress;
                tDisp(istre,:)      = vars.d_u;
                obj.assignVarsToPrint(istre);
            end
%             nstre = obj.material.nstre;
%             basis = diag(ones(nstre,1));
%             Ch = zeros(nstre,nstre);
%                     nelem = size(obj.material.C,3);
%                     dim = obj.getDimensions();
%                     npnod = dim.nnodes;
%                     ndofs = dim.ndimf;
%                     ngaus = obj.quadrature.ngaus;
%                     tStrn  = zeros(nstre,ngaus,nstre,nelem);
%                     tStrss = zeros(nstre,ngaus,nstre,nelem);
%                     tDisp  = zeros(nstre,ndofs);
%                     for istre=1:nstre
%                         obj.vstrain = basis(istre,:);
%                         obj.solve();
%                         vars = obj.computeStressStrainAndCh();
%                         Ch(:,istre)         = vars.stress_homog;
%                         tStrn(istre,:,:,:)  = vars.strain;
%                         tStrss(istre,:,:,:) = vars.stress;
%                         tDisp(istre,:)      = vars.d_u;
%                       %  obj.assignVarsToPrint(istre);
%                     end
                    obj.variables.Chomog  = Ch;
                    obj.variables.tstrain = tStrn;
                    obj.variables.tstress = tStrss;
                    obj.variables.tdisp   = tDisp;
        end

    end

    methods (Access = private)
        function vars = computeStressStrainAndCh(obj)
            vStrn       = obj.vstrain;
            vars        = obj.variables;
            Cmat        = obj.material.C;
            nstre       = obj.material.nstre;
            nelem       = size(Cmat,3);
            ngaus       = obj.quadrature.ngaus;
            dV          = obj.mesh.computeDvolume(obj.quadrature)';
            strainFluct = vars.strain;
            stressFluct = vars.stress;
            
            stress      = zeros(ngaus,nstre,nelem);
            strain      = zeros(ngaus,nstre,nelem);
            stressHomog = zeros(nstre,1);
            
            for igaus = 1:ngaus
                strain(igaus,1:nstre,:) = vStrn.*ones(1,nstre,nelem) + strainFluct(igaus,1:nstre,:);
                for istre = 1:nstre
                    for jstre = 1:nstre
                        Cij  = squeeze(Cmat(istre,jstre,:,igaus));
                        C    = squeeze(Cij);
                        strs = squeeze(stress(igaus,istre,:));
                        strn = squeeze(strain(igaus,jstre,:));
                        stress(igaus,istre,:) = strs + C.* strn;
                    end
                    strs               = squeeze(stress(igaus,istre,:));
                    stressHomog(istre) = stressHomog(istre) + (strs)'*dV(:,igaus);
                end
            end

            vars.stress_fluct = stressFluct;
            vars.strain_fluct = strainFluct;

            vars.stress       = stress;
            vars.strain       = strain;
            vars.stress_homog = stressHomog;
            obj.variables     = vars;
        end

        function assignVarsToPrint(obj, istre)
            vars = obj.variables;
            ndimField = obj.displacementFun.ndimf; 
            obj.variables2print{istre}.d_u          = reshape(vars.d_u',ndimField,[])';
            obj.variables2print{istre}.fext         = reshape(vars.fext',ndimField,[])';
            obj.variables2print{istre}.stress       = vars.stress;
            obj.variables2print{istre}.strain       = vars.strain;
            obj.variables2print{istre}.stress_fluct = vars.stress_fluct;
            obj.variables2print{istre}.strain_fluct = vars.strain_fluct;
        end
    end
end