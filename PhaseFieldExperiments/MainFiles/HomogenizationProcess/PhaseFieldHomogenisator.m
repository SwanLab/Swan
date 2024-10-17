classdef PhaseFieldHomogenisator < handle

    properties (Access = private)
        E
        nu
    end

    methods (Access = public)
        
        function obj = PhaseFieldHomogenisator(obj)
            cParams.E = 210;
            cParams.nu = 0.3;
            obj.init(cParams);
        end
        
        function [mat,phi] = computeHomogMaterial(obj,holeType,aType,steps)
            holeMax = obj.computeHoleMax(holeType);
            holeLength = linspace(1e-10,holeMax,steps);
            mat = zeros(3,3,length(holeLength));
            phi = zeros(length(holeLength),1);
            for i=1:steps
                l = holeLength(i);
                mat(:,:,i) = obj.computeHomogenisation(l,holeType);
                phi(i) = obj.computeDissipationMetric(l,holeType,aType);
            end
            %Chomog = obj.createMaterialFun(mat,alpha);
        end
        
        function [mat,holeLength] = computeIsotropicMaterial(obj,aType,steps)
            constant =obj.E/(1-obj.nu^2);
            C(1,1) = constant;
            C(1,2) = constant*obj.nu;
            C(2,1) = constant*obj.nu;
            C(2,2) = constant;
            C(3,3) = constant*(1-obj.nu)/2;
            
            phi = linspace(0,1,steps);
            switch aType
                case "AT1"
                    alpha = phi;
                case "AT2"
                    alpha = phi.^2;
            end
            Ciso = cell(3,3);
            for i=1:3
                for j=1:3
                    sM.coord = alpha';
                    sM.connec = (1:length(alpha)-1)' + [0,1];
                    s.mesh = Mesh.create(sM);
                    s.fValues = ((1-phi).^2*C(i,j))';
                    s.order = 'P1';
                    Ciso{i,j} = LagrangianFunction(s);
                end
            end
            holeLength = linspace(0,1,steps);
            mat = zeros(3,3,length(phi));
            mat(1,1,:) = constant*(1-phi).^2 ;
            mat(1,2,:) = constant*obj.nu*(1-phi).^2;
            mat(2,1,:) = constant*obj.nu*(1-phi).^2;
            mat(2,2,:) = constant*(1-phi).^2;
            mat(3,3,:) = constant*(1-obj.nu)*(1-phi).^2 /2;
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.E = cParams.E;
            obj.nu = cParams.nu;
        end
        
        function holeMax = computeHoleMax(obj,holeType)
            switch holeType
                case "Circle"
                    holeMax = 0.5;
                case "Square"
                    holeMax = 0.99;
                case "Full"
                    holeMax = 1;
            end
        end
        
        function matHomog = computeHomogenisation(obj,l,holeType)
            mat = Tutorial02p2FEMElasticityMicro(l,holeType,obj.E,obj.nu);
            matHomog = mat.stateProblem.Chomog;
            
            figure(1)
            cla reset
            mat.mesh.plot
        end
        
        function alpha = computeDissipationMetric(obj,l,holeType,aType)
            switch aType
                case "Area"
                    switch holeType
                        case "Circle"
                            alpha = (pi*l^2)/(pi*0.5*0.5);
                        case "Square"
                            alpha = l^2;
                        case "Full"
                            alpha = 0;
                    end
                case "Perimeter"
                    switch holeType
                        case "Circle"
                            alpha = 2*pi*l/pi;
                        case "Square"
                            alpha = 4*l/4;
                        case "Full"
                            alpha = 0;
                    end
            end
        end
        
        function Chomog = createMaterialFun(obj,mat,alpha)
            Chomog = cell(3,3);
            for i=1:3
                for j=1:3
                    sM.coord = alpha;
                    sM.connec = [1:length(alpha)-1]' + [0,1];
                    s.mesh = Mesh.create(sM);
                    s.fValues = squeeze(mat(i,j,:));
                    s.order = 'P1';
                    Chomog{i,j} = LagrangianFunction(s).project('P2');
                end
            end
        end
        
    end
    
end