classdef SystemSolver < handle

    properties (Access = public)
        displacements
        force
        stiffnessMatrix
        volume
        charFunc
    end

    properties (Access = private)
        mesh
        pdeCoeff
        tensor
        bc
        E
        gamma
        levelSet
        effectiveTensor
    end

    methods (Access = public)

        function obj = SystemSolver(cParams)
            obj.init(cParams)
            obj.computeGamma();
            obj.computeCharacteristicFunction();
            obj.computeEffectiveTensor();
        end

        function [U,F,vol] = computeStiffnessMatrixAndForce(obj)
            p = obj.mesh.p;
            t = obj.mesh.t;
            a = obj.pdeCoeff.a;
            f = obj.pdeCoeff.f;
            c = obj.effectiveTensor;
            gamma = obj.gamma;
            tfi = obj.charFunc;
    
          %  tgamma = pdeintrp(p,t,fi*gamma');  %for P1-projection aproach -- WHAT IS THIS?
            [K,~,F] = assema(p,t,c,a,f);
            [K,F] = pdeupdate(K,F,obj.bc,obj.mesh);
            U = K \ F;  % solve linear system
            %     vol  = area*tchi';  %calculate volume assigned to each material including void -- P1 projection aproach
            vol  = obj.mesh.area*tfi';  %calculate volume assigned to each material including void -- mixed formulation aproach


            % assembly = assema(s); % s'ha de construir la classe!!
            %  obj.stiffnessMatrix = assembly.computeStiffMatrix();
            %  obj.force = assembly.computeForce();


        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.E(1) = cParams.matProp.matA.young;
            obj.E(2) = cParams.matProp.matB.young;
            obj.E(3) = cParams.matProp.matC.young;
            obj.E(4) = cParams.matProp.matD.young;

            obj.mesh     = cParams.mesh;
            obj.pdeCoeff = cParams.pdeCoeff;
            obj.bc       = cParams.bc;
            obj.levelSet = cParams.psi;
        end
        
        function computeGamma(obj)
            obj.gamma = obj.E./obj.E(1); % contrast for each material
        end


        function computeCharacteristicFunction(obj)
            s.p = obj.mesh.p;
            s.t = obj.mesh.t;
            s.pdeCoeff = obj.pdeCoeff; 
            s.psi = obj.levelSet;

            chi = CharacteristicFunctionComputer(s); % s'ha de construir la classe - charfunc!!
          %  [ ~,tfi ] = charfunc(p,t,psi);
            obj.charFunc = chi.tfi;
        end

        function computeEffectiveTensor(obj)
            tgamma = obj.gamma*obj.charFunc; % for mixed formulation approach
            c0 =  obj.pdeCoeff.tensor(:,1);
            obj.effectiveTensor = c0*tgamma; % effective elasticity tensor
        end

    


    end
end
   
%

