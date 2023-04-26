classdef ComplianceRobustComputer < handle
    properties (Access = public)
        E
        I
        D
    end
    properties (Access = private)
        iterations
        mesh
        structure
        projectorParameters
        filterParameters
        solverParameters
        data
    end
    methods (Access = public)
        function obj = ComplianceRobustComputer(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.computeInitialParameters();
            obj.computeCost();
            obj.optimize();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.iterations = cParams.iterations;
            obj.readDataFile(cParams.meshName);
        end 
        function readDataFile(obj,meshFileName)
            run(meshFileName);
            obj.data.Data_prb = Data_prb;
            obj.data.External_border_elements = External_border_elements;
            obj.data.External_border_nodes = External_border_nodes;
            obj.data.gidcoord = gidcoord;
            obj.data.gidlnods = gidlnods;
            obj.data.lnodes = lnodes;
            obj.data.meshFileName = meshFileName;
            obj.data.nodesolid = nodesolid;
            obj.data.pointload_complete = pointload_complete;           
        end
        function computeInitialParameters(obj)
            obj.computeMeshParameters();
            obj.computeStructureParameters();
            obj.computeSolverParameters();
            obj.computeProjectorParameters
            obj.computeFilterParameters();
            obj.computeInitialFields();
            obj.computeInitialVolumens();
        end
        function computeCost(obj)          
            %Get intial cost
            s.mesh = obj.mesh; 
            s.structure = obj.structure;
            
            s.designField = obj.E.designField;
            obj.E.designCost = DesignCost(s);
            obj.E.designCost.computeCost();
            s.designField = obj.I.designField;
            obj.I.designCost = DesignCost(s);
            obj.I.designCost.computeCost()
            s.designField = obj.D.designField;
            obj.D.designCost = DesignCost(s);
            obj.D.designCost.computeCost()


        end 
        function optimize(obj)
            %Optimización
            s.mesh = obj.mesh;
            s.structure=obj.structure;
            s.structure.elementType = 'TRIANGULAR';
            s.projectorParameters = obj.projectorParameters;
            s.filterParameters = obj.filterParameters;
            s.solverParameters =obj.solverParameters;
            s.solverParameters.initialCost = obj.E.designCost.cost;
            s.iterations = obj.iterations;
            s.E = obj.E;
            s.I = obj.I;
            s.D = obj.D;
            B = Optimizer(s);
            B.compute();
            obj.E = B.E;
            obj.I = B.I;
            obj.D = B.D;
        end 
        function computeMeshParameters(obj)
            %Mesh parameters
            %Malla estandar 400x80
            obj.mesh.neumanCondition       = -1e-3;
            obj.mesh.output    = 2;
            %Compute Elemental stifness Matrix
           

%           obj.mesh.elementNumberX    = length(unique(obj.data.gidcoord(:,2)));
%           obj.mesh.elementNumberY    = length(unique(obj.data.gidcoord(:,3)));
%           obj.mesh.elementNumberX    = 160;
%           obj.mesh.elementNumberY    = 80;
            s.elementType = 'TRIANGULAR';
            s.data = obj.data;
%             s.elementNumberX =  obj.mesh.elementNumberX;
%             s.elementNumberY =  obj.mesh.elementNumberY;
            B = GeometryComputer(s);
            B.compute();
            obj.mesh.degress.all   = B.degress.all;
            obj.mesh.degress.free  = B.degress.free;
            obj.mesh.degress.fixed = B.degress.fixed;
            obj.mesh.conectivityMatrixMat = B.conectivityMatrixMat;
        end
        function computeStructureParameters(obj)
             obj.structure.t       = 1;
             %Parámetros del material:
            obj.structure.elasticModuleNeutral      = 1;
            obj.structure.elasticModuleMinimun    = 1e-6;
            obj.structure.poissonCoefficient      = 0.3;
            obj.structure.penalization   = 3;

            %Compute Elemental stifness Matrix
            s.elementType = 'SQUARE';
            %s.elementType = cell2mat(obj.Data.Data_prb(1));
            s.t = obj.structure.t;
            s.poissonCoefficient =obj.structure. poissonCoefficient;
            B = ElementalStiffnessMatricesComputer(s);
            B.compute();
            obj.structure.elementalStiffnessMatrix = B.elementalStiffnessMatrix;
        end 
        function computeSolverParameters(obj)
            obj.solverParameters.costChange  = 1;
            obj.solverParameters.minDensity    = zeros(obj.mesh.elementNumberY,obj.mesh.elementNumberX);
            obj.solverParameters.maxDensity    =  ones(obj.mesh.elementNumberY,obj.mesh.elementNumberX);
            obj.solverParameters.minDensity =  obj.solverParameters.minDensity(:);
            obj.solverParameters.maxDensity =  obj.solverParameters.maxDensity(:);
            
            obj.solverParameters.numberRestriction       = 4;
            obj.solverParameters.variableNumber       = obj.mesh.elementNumberX*obj.mesh.elementNumberY;
            obj.solverParameters.low     = obj.solverParameters.minDensity;
            obj.solverParameters.upp     = obj.solverParameters.maxDensity;
            obj.solverParameters.a0      = 1;
            obj.solverParameters.mmaParameter.a       = [1 1 1 0]';
            obj.solverParameters.mmaParameter.e       = 1000*ones(obj.solverParameters.numberRestriction,1);
            obj.solverParameters.mmaParameter.d       = 0*ones(obj.solverParameters.numberRestriction,1);
            obj.solverParameters.xold1   = [];
            obj.solverParameters.xold2   = [];
        end 
        function computeFilterParameters(obj)
            %Filter parameters:
            obj.filterParameters.volumenFraction =   0.5;
            obj.filterParameters.minimunInfluenceRadios    =   1.5;
            
            s.elementNumberX =  obj.mesh.elementNumberX;
            s.elementNumberY =  obj.mesh.elementNumberY;
            s.minimunInfluenceRadios = obj.filterParameters.minimunInfluenceRadios;
            B = weightFilterComputer(s);
            B.compute();
            obj.filterParameters.H = B.H;
            obj.filterParameters.Hs = B.Hs;
        end

        function computeProjectorParameters(obj)
            eta      = 0.25;
            obj.projectorParameters.eta.E     = 1-eta;
            obj.projectorParameters.eta.I     = 0.5;
            obj.projectorParameters.eta.D     = eta;
            obj.projectorParameters.beta     = 1;
        end
        function computeInitialVolumens(obj)
            s.mesh = obj.mesh;
            s.filterParameters =obj.filterParameters;            
            s.designField =obj.E.designField;
            obj.E.designVolumen = DesignVolumen(s);
            obj.E.designVolumen.computeVolumenFraction(obj.D,obj.I);
            s.designField =obj.I.designField;            
            obj.I.designVolumen = DesignVolumen(s);
            obj.I.designVolumen.computeVolumenFraction(obj.D,obj.I);
            s.designField =obj.D.designField;         
            obj.D.designVolumen = DesignVolumen(s);
            obj.D.designVolumen.computeVolumenFraction(obj.D,obj.I);


        end 

        function computeInitialFields(obj)
            s.field = obj.filterParameters.volumenFraction*ones(obj.mesh.elementNumberY,obj.mesh.elementNumberX);
            s.field([],[]) = 1;
            s.filterParameters =  obj.filterParameters;
            s.projectorParameters.beta =  obj.projectorParameters.beta;
            s.mesh = obj.mesh;

            s.projectorParameters.eta =  obj.projectorParameters.eta.E; 
            obj.E.designField = DesignField(s);
            s.projectorParameters.eta =  obj.projectorParameters.eta.I;
            obj.I.designField =  DesignField(s);
            s.projectorParameters.eta =  obj.projectorParameters.eta.D;
            obj.D.designField =  DesignField(s);

            obj.E.designField.filteredField = obj.E.designField.field;
            obj.I.designField.filteredField =  obj.I.designField.field;
            obj.D.designField.filteredField = obj.D.designField.field;

            obj.E.designField.project();
            obj.I.designField.project();
            obj.D.designField.project();


        end
        
    end
end