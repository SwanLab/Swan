classdef NewElasticProblemMicro < handle
    
    properties (Access = public)
        variables
        uFun, strainFun, stressFun
    end

    properties (Access = private)
        elasticParams
        material
        quadrature
    end

    methods (Access = public)

        function obj = NewElasticProblemMicro(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function v = computeGeometricalVolume(obj)
            v = 1;%sum(sum(obj.geometry.dvolu));
        end

        function dvolu = getDvolume(obj)
            dvolu  = obj.elasticParams.mesh.computeDvolume(obj.quadrature);
        end

        function setMatProps(obj,s)
           obj.material.compute(s);
        end

        function mesh = getMesh(obj)
            mesh  = obj.elasticParams.mesh;
        end
        
        function interp = getInterpolation(obj)
            interp  = obj.elasticParams.mesh.interpolation;
            interp.computeShapeDeriv(obj.quadrature.posgp);
        end

        function quad = getQuadrature(obj)
            quad = obj.quadrature;
        end

        function setC(obj, C)
            obj.material.C = C;
        end

        function dim = getDimensions(obj)
            strdim = regexp(obj.elasticParams.dim,'\d*','Match');
            nDimf  = str2double(strdim);
            d.ndimf  = nDimf;
            d.nnodes = obj.elasticParams.mesh.nnodes;
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.elasticParams.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function obj = computeChomog(obj)
            nVoigt = obj.material.nstre;
            nGaus = obj.quadrature.ngaus;
            nElem = size(obj.material.C,3);
            nPnod = obj.elasticParams.mesh.nnodes; %size(obj.displacementFun.fValues,1);
            strdim = regexp(obj.elasticParams.dim,'\d*','Match');
            nDimf  = str2double(strdim);
            nDofs = nPnod*nDimf;
            
            tStrn  = zeros(nVoigt,nGaus,nVoigt,nElem);
            tStrss = zeros(nVoigt,nGaus,nVoigt,nElem);
            tDisp  = zeros(nVoigt,nDofs);
            Ch = zeros(nVoigt,nVoigt);
            
            params  = obj.elasticParams;
            prob = NewElasticProblem(params);

            for iVoigt = 1:nVoigt
                vstrain = obj.computeVstrain(iVoigt);
                prob = obj.solveElasticProblem(prob,vstrain);
                vars = obj.computeStressStrainFluctuations(prob, vstrain);
                Ch(:,iVoigt)         = vars.stress_homog;
                tStrn(iVoigt,:,:,:)  = vars.strain;
                tStrss(iVoigt,:,:,:) = vars.stress;
                tDisp(iVoigt,:)      = prob.variables.d_u;
                obj.uFun{iVoigt}      = prob.uFun;
                obj.strainFun{iVoigt} = prob.strainFun;
                obj.stressFun{iVoigt} = prob.stressFun;
            end
            obj.variables.Chomog  = Ch;
            obj.variables.tstrain = tStrn;
            obj.variables.tstress = tStrss;
            obj.variables.tdisp   = tDisp;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.elasticParams = cParams;
            obj.material = cParams.material;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.elasticParams.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function vstrain = computeVstrain(obj, iVoigt)
            nVoigt  = obj.material.nstre;
            basis   = diag(ones(nVoigt,1));
            vstrain = basis(iVoigt,:);
        end

        function prob = solveElasticProblem(obj, prob, vstrain)
            prob.setVstrain(vstrain);
            prob.solve();
        end

        function vars = computeStressStrainFluctuations(obj, prob, vstrain)
%             vars  = obj.variables;
            Cmat  = obj.material.C;
            nstre = obj.material.nstre;
            nelem = size(Cmat,3);
            ngaus = obj.quadrature.ngaus;
            dV = obj.elasticParams.mesh.computeDvolume(obj.quadrature)';
            strainFluct = permute(prob.strainFun.fValues, [2 1 3]);
            stressFluct = permute(prob.stressFun.fValues, [2 1 3]);
            
            stress = zeros(ngaus,nstre,nelem);
            strain = zeros(ngaus,nstre,nelem);
            stressHomog = zeros(nstre,1);
            
            for igaus = 1:ngaus
                strain(igaus,1:nstre,:) = vstrain.*ones(1,nstre,nelem) + strainFluct(igaus,1:nstre,:);
                for istre = 1:nstre
                    for jstre = 1:nstre
                        Cij  = squeeze(Cmat(istre,jstre,:,igaus));
                        C    = squeeze(Cij);
                        strs = squeeze(stress(igaus,istre,:));
                        strn = squeeze(strain(igaus,jstre,:));
                        stress(igaus,istre,:) = strs + C.* strn;
                    end
                    strs = squeeze(stress(igaus,istre,:));
                    stressHomog(istre) = stressHomog(istre) + (strs)'*dV(:,igaus);
                end
            end

            vars.stress_fluct = stressFluct;
            vars.strain_fluct = strainFluct;

            vars.stress = stress;
            vars.strain = strain;
            vars.stress_homog = stressHomog;
            obj.variables = vars;
        end
    end

end
