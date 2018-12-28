classdef testFEMPrinting < ...
            testNotShowingError & testFemComputation
    
    properties (Access = protected)
        testName = 'test2d_quad';  
        postProcessor = 'Elasticity';
        filesHaveChanged
        fileOutputName
    end
    
    methods (Access = public)
        
        function obj = testFEMPrinting()
            obj.init()
            obj.print()
            obj.compareFiles()
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            hasPassed = ~obj.filesHaveChanged;
        end
        
        function selectComputedVar(obj)
        end
        
    end

    methods (Access = private)
        
        function init(obj)
            obj.fileOutputName = ['testFemPrinting'];
        end
        
        function compareFiles(obj)
            hasMshChanged = obj.compareFile('.msh');
            hasResChanged = obj.compareFile('.res');
            obj.filesHaveChanged = hasMshChanged || hasResChanged;
        end
        
        function hasChanged = compareFile(obj,extension)
            out   = obj.fileOutputName;
            name1 = ['tests/PrintingTests/PrintedFiles/',out,'_u_1.flavia',extension];
            name2 = ['Output/',out,'/',out,'1.flavia',extension];
            command = ['diff ', name1, ' ', name2];
            [hasChanged,~] = system(command);
        end        
        
        
        function print(obj)            
            postprocess = Postprocess(obj.postProcessor);
            d = obj.createPostProcessDataBaseStructre();           
            postprocess.print(d);
        end
        
        function d = createPostProcessDataBaseStructre(obj)
            d.fields   = obj.fem.variables;
            d.iter     = 0;
            d.outFileName  = obj.fileOutputName;            
            d.nfields = 1;
            d.coordinates = obj.fem.element.interpolation_u.xpoints;
            d.connectivities = obj.fem.element.interpolation_u.T;
            d.ngaus = obj.fem.element(1).quadrature.ngaus;
            d.posgp = obj.fem.element(1).quadrature.posgp';
            d.nnode = obj.fem.element.nnode;
            d.npnod = obj.fem.element.interpolation_u.npnod;
            d.gtype = obj.fem.mesh.geometryType;
            d.ndim  = obj.fem.element.interpolation_u.ndime;
            d.pdim  = obj.fem.mesh.pdim;
            d.ngaus = obj.fem.element(1).quadrature.ngaus;
            d.posgp = obj.fem.element(1).quadrature.posgp';
            d.ptype = obj.fem.mesh.ptype;
            switch  d.gtype
                case 'TRIANGLE'
                    d.etype = 'Triangle';
                case 'QUAD'
                    d.etype = 'Quadrilateral';
                case 'TETRAHEDRA'
                    d.etype = 'Tetrahedra';
                case 'HEXAHEDRA'
                    d.etype = 'Hexahedra';
            end
            d.nelem    = obj.fem.element.nelem;
        end
    end
    
end

