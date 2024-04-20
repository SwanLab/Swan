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
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.E(1) = cParams.A.young;
            obj.E(2) = cParams.B.young;
            obj.E(3) = cParams.C.young;
            obj.E(4) = cParams.D.young;

            obj.mesh = cParams.mesh;
            obj.pdeCoeff = cParams.pdeCoeff;

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

            tfi = CharacteristicFunctionComputer(s); % s'ha de construir la classe - charfunc!!
            obj.charFunc = tfi.compute();
        end

        function computeEffectiveTensor(obj)
            tgamma = obj.gamma*obj.charFunc; % for mixed formulation approach
            c0 =  obj.pdeCoeff.tensor(:,1);
            obj.effectiveTensor = c0*tgamma; % effective elasticity tensor
        end

        function computeStiffnessMatrixAndForce(obj)
            s.p = obj.mesh.p;
            s.t = obj.mesh.t;
            s.a = obj.pdeCoeff.a;
            s.f = obj.pdeCoeff.f;
            s.c = obj.pdeCoeff.tensor; 

            assembly = assema(s); % s'ha de construir la classe!!
            obj.stiffnessMatrix = assembly.computeStiffMatrix();
            obj.force = assembly.computeForce();

        end


    end
end
   
% tgamma = pdeintrp(p,t,fi*gamma');  %for P1-projection aproach -- WHAT IS THIS? 
    [K,~,F] = assema(p,t,c,a,f);
    [K,F] = pdeupdate(K,F,bc,mesh);
    U = K \ F;  % solve linear system
%     vol  = area*tchi';  %calculate volume assigned to each material including void -- P1 projection aproach
    vol  = area*tfi';  %calculate volume assigned to each material including void -- mixed formulation aproach


