classdef testTopOptPrinting < testNotShowingError ...
                                  & testTopOptComputation ...
                                  & testPrintingDescriptor
    
    properties (Access = protected)
      filesHaveChanged
      iter
    end
    
    properties (Access = protected, Abstract)
       postProcessor 
    end
    
    methods (Access = protected)
        
        function obj = testTopOptPrinting()
            obj.init();
            obj.print();
            obj.compareFiles();
        end
        
        function init(obj)
            obj.iter = 0;
        end
        
        function print(obj)
            mesh      = obj.topOpt.mesh;
            d.fields  = obj.topOpt.x;
            d.outFileName = obj.fileOutputName;
            d.iter    = obj.iter;
            
            d.coordinates = mesh.coord;
            d.connectivities = mesh.connec;
            d.nnode = size(mesh.connec,2);
            d.npnod = size(mesh.coord,1);  % Number of nodes

            d.gtype = mesh.geometryType;
            d.pdim  = mesh.pdim;
            switch d.pdim
                case '2D'
                    d.ndim=2;
                case '3D'
                    d.ndim=3;
            end
            d.ptype = mesh.ptype;
            
            switch  d.gtype %gid type
                case 'TRIANGLE'
                    d.etype = 'Triangle';
                case 'QUAD'
                    d.etype = 'Quadrilateral';
                case 'TETRAHEDRA'
                    d.etype = 'Tetrahedra';
                case 'HEXAHEDRA'
                    d.etype = 'Hexahedra';
            end
            d.nelem = size(mesh.connec,1); % Number of elements            
            postprocess = Postprocess(obj.postProcessor);
            postprocess.print(d);
        end
        
    end
    
end

