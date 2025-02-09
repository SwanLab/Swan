classdef DehomogenisationPrinter < handle
    
    properties (Access = private)
        mesh
        dataRes
        alpha
        density
    end
    
    properties (Access = private)
        fileName
        folderPath        
    end
    
    methods (Access = public)
        
        function obj = DehomogenisationPrinter(cParams)
            obj.init(cParams);
        end
        
        function print(obj)
            obj.wrapResMeshData();
            obj.projectAlphaToNodes();
            obj.projectDensityToNodes();
            obj.printDehomogResFile();
            obj.printDehomogMshFile();            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName   = cParams.fileName;
            obj.folderPath = cParams.folderPath;
        end

        function wrapResMeshData(obj)
            s.fileName   = obj.fileName;
            s.folderPath = obj.folderPath;
            w = WrapperMshResFiles(s);
            w.compute();
            obj.dataRes = w.dataRes;
            obj.mesh    = w.mesh;            
        end
        
        function projectAlphaToNodes(obj)
            alphaN = obj.dataRes.AlphaGauss;
            alpha1 = obj.projectInNodesScalarPieceWiseFunction(alphaN(:,1));
            alpha2 = obj.projectInNodesScalarPieceWiseFunction(alphaN(:,2));
            normA = sqrt(alpha1.^2 + alpha2.^2);
            obj.alpha = [alpha1,alpha2]./normA;
        end

        function projectDensityToNodes(obj)
            alphaN = obj.dataRes.DensityGauss;
            dens   = obj.projectInNodesScalarPieceWiseFunction(alphaN(:,1));
            obj.density = dens;
        end
        
        function alpha = projectInNodesScalarPieceWiseFunction(obj,fValues)
            s.mesh    = obj.mesh;
            s.fValues = fValues;
            s.order   = 'P0';
            f = P0FunctiLagrangianFunctionon(s);
            alpha = f.projectToLinearNodalFunction();
        end
        
        function printDehomogResFile(obj)
            fD = [obj.fileName,'DehomogRes.txt'];
            fOutName = fullfile(obj.folderPath,fD);
            s.fileName = fOutName;
            s.values(:,1) = obj.dataRes.DesignVar1;
            s.values(:,2) = obj.dataRes.DesignVar2;
            s.values(:,3) = obj.alpha(:,1);
            s.values(:,4) = obj.alpha(:,2);
            s.values(:,5) = obj.dataRes.SuperEllipseExponent;
            s.values(:,6) = obj.density;            
            p = DehomogenisationOutputPrinter(s);
            p.print();
        end
        
        function printDehomogMshFile(obj)
            fD = [obj.fileName,'DehomogMesh.txt'];
            fOutName = fullfile(obj.folderPath,fD);
            s.fileName = fOutName;
            s.mesh     = obj.mesh;
            p = DehomogenisationMeshPrinter(s);
            p.print();
        end
        
    end
    
end