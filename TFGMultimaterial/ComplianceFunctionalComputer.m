classdef ComplianceFunctionalComputer < handle

    properties (Access = public)
        J
        dJ
    end

    properties (Access = private)
        compliance
        materialInterpolator
        material
        designVariable
        F
        U
        energy0
        nMat
        tgamma
        mesh
        bc
        pdeCoeff
        mat
        young
        gamma
        tE
        beta
        alpha
        la0
        mu0
        UFun
    end

    methods (Access = public)
        
        function obj = ComplianceFunctionalComputer(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            obj.material.setDesignVariable(x);
            C = obj.material.obtainTensor();
            obj.solveFEM(C);
            obj.computeGamma();
            obj.computeTgamma(x);

            for i = 1:obj.nMat
                for j = 1:obj.nMat
                    obj.materialInterpolator.computeFromTo(i,j);
                    dCfun = obj.material.obtainTensorDerivative();
                    dCeval =  squeezeParticular(dCfun.evaluate([0;0]),3);
                    dC(:,:,:,i,j) = dCeval;
                    dC(1,2,:,i,j) = 0;
                    dC(2,1,:,i,j) = 0;
                    dC(2,3,:,i,j) = 0;
                    dC(3,2,:,i,j) = 0;
                    dC(1,3,:,i,j) = dCeval(1,2,:);
                    dC(2,2,:,i,j) = 4*dCeval(3,3,:);
                    dC(3,1,:,i,j) = dCeval(2,1,:);
                    dC(3,3,:,i,j) = dCeval(2,2,:);
                end
            end

            J = obj.computeFunction();
            dJ = obj.computeGradient(dC,x);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.compliance = cParams.complianceFromConstitutive;
            obj.materialInterpolator = cParams.materialInterpolator;
            obj.material = cParams.material;
            obj.nMat = cParams.nMat;
            obj.mat = cParams.mat;
            obj.mesh = cParams.mesh;
            obj.pdeCoeff = cParams.pdeCoeff;
            obj.bc = cParams.bc;
        end

        function J = computeFunction(obj)
            J = 0.5*dot(obj.F,obj.U); 
            if isempty(obj.energy0)
                obj.energy0 = J;
            end
            J = J/obj.energy0;
        end

        function dJ = computeGradient(obj,dC,x)
            e = obj.computeStrain();

            for i = 1:obj.nMat
                for j = 1:obj.nMat
                    % for z=1:size(dC,3)
                    %     derTop(z) = e(:,z)'*dC(:,:,z,i,j)*e(:,z);
                    % end
                    
                    dCij = dC(:,:,:,i,j);
                    ev(:,1,:) = e;
                    s = pagemtimes(dCij,ev);
                    derTop = squeezeParticular(sum(ev.*s,1),1);
                    %derTop2 = e(:,z)'*dC(:,:,z,i,j)*e(:,z)
                    
                    
                    
                    TD{i,j} = derTop; % lo que abans era TD
                    DJ{i,j} = TD{i,j}/obj.energy0; % adimensionat
                end
            end

            dJ = obj.smoothGradient(DJ,x);

            dJ = reshape(dJ,[],1);
        end

        function C = computeElasticTensor(obj,x)
            s.matProp           = obj.mat;
            s.pdeCoeff          = obj.pdeCoeff;
            s.bc                = obj.bc;
            s.mesh                 = obj.mesh;
            s.designVariable    = x;
            
            constituitiveTensor = ElasticTensorComputer(s);
            C = constituitiveTensor.C;
        end

        function solveFEM(obj,C)
            s.mesh               = obj.mesh;
            s.scale              = 'MACRO';
            s.material           = C;
            s.dim                = '2D';
            s.boundaryConditions = obj.bc;
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';
            s.solverCase         = 'DIRECT';

            fem         = ElasticProblem(s);
            fem.solve();
            displ    = fem.uFun.fValues; 
            %obj.displFun = fem.uFun;
            obj.U          = reshape(displ, [size(displ,1)*2, 1]);
            force      = fem.forces; 
            Fx         = force(1:2:end); % Fx values are at odd indices
            Fy         = force(2:2:end); % Fy values are at even indices
            forcesVect = [Fx Fy];
            obj.F          = reshape(forcesVect, [size(displ,1)*2, 1]);
            obj.UFun = fem.uFun;
        end

        function computeGamma(obj)
            E1 = obj.mat.A.young;
            E2 = obj.mat.B.young;
            E3 = obj.mat.C.young;
            E4 = obj.mat.D.young;

            obj.young = [E1 E2 E3 E4];

            obj.gamma(1) = E1/E1;
            obj.gamma(2) = E2/E1;
            obj.gamma(3) = E3/E1;
            obj.gamma(4) = E4/E1;    
        end

        function computeTgamma(obj,x)
            charfun = x.obtainDomainFunction(); 
            [~,tfi] = charfun.computeAtNodesAndElements();

            obj.tgamma = obj.gamma*tfi.fValues'; %Mixed formulation method
            E1 = obj.mat.A.young;
            obj.tE = E1*obj.tgamma; 
        end

        function computeBetaAndAlpha(obj)
            nu = obj.mat.A.nu;
            obj.beta = (1+nu)/(1-nu); 
            obj.alpha = (3-nu)/(1+nu);
        end

        function computeLameParameters(obj)
            nu = obj.mat.A.nu;
            E1 = obj.mat.B.young;
        
            obj.la0 = nu.*E1./((1+nu).*(1-2.*nu));
            obj.mu0 = E1./(2.*(1+nu)); % plane strain
            obj.la0 = 2.*obj.mu0.*obj.la0./(obj.la0+2.*obj.mu0); % plane stress 

        end

        function e =computeStrain(obj)
            tgamma3 = [obj.tgamma;obj.tgamma;obj.tgamma]; 
            uFun = obj.UFun;
            str = SymGrad(uFun);
            strainSw = str.evaluate([0;0]);
            strain(1,:) = strainSw(1,:);
            strain(2,:) = strainSw(3,:)/2;
            strain(3,:) = strainSw(2,:);
            e = strain.*tgamma3;
        end
        
        function dJ = smoothGradient(obj,TD,x)
            s.mesh           = obj.mesh;
            s.designVariable = x;
            multGrad         = MultimaterialGradientComputer(s);
            dt               = multGrad.compute(TD);
            t = obj.mesh.connec';
            p = obj.mesh.coord';
            dJ = pdeprtni(p,t,dt);
        end

     

        
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end


end