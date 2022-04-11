classdef MatlabFileWriter < handle
    
    properties (Access = public)
        filename
        coord
        connec
        nodes
        masterSlaveIndex
    end
    
    methods (Access = public)
        
        function obj = MatlabFileWriter(cParams)
            obj.init(cParams);
        end
        
        function write(obj)
            obj.writeMeshFile();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filename = cParams.filename;
            obj.coord = cParams.coord;
            obj.connec = cParams.connec;
            obj.nodes = cParams.nodes;
            obj.masterSlaveIndex = cParams.masterSlaveIndex;
        end
        
        function writeMeshFile(obj)
            Data_prb = {'''TRIANGLE''','''SI''','''2D''','''Plane_Stress''','''ELASTIC''','''MICRO'''};
            % Initialization of the document
            fileID = fopen(obj.filename,'w');
            fprintf(fileID,'%c','%% Data file');
            fprintf(fileID,'\n\n');
            
            % Characteristics
            fprintf(fileID,'%c','Data_prb = {');
            fprintf(fileID,'\n');
            for k = 1:length(Data_prb)
                fprintf(fileID,'%s ;\r\n',Data_prb{k});
            end
            fprintf(fileID,'%c','};');
            fprintf(fileID,'\n\n\n');
            
            % Coordinates
            fprintf(fileID,'%c','coord = [');
            fprintf(fileID,'\n');
            for k = 1:size(obj.coord,1)
                fprintf(fileID,'%g %g %g %g\r\n',k,obj.coord(k,1),obj.coord(k,2),0);
            end
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n');
            
            % Point loads
            fprintf(fileID,'%c','pointload = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n\n');
            
            % Connectivities
            fprintf(fileID,'%c','connec = [');
            fprintf(fileID,'\n');
            for k = 1:size(obj.connec,1)
                fprintf(fileID,'%g %g %g %g\r\n',k,obj.connec(k,1),obj.connec(k,2),obj.connec(k,3));
            end
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n');
            
            % Variable prescribed
            fprintf(fileID,'%c','%% Variable Prescribed');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s \t %s \t %s','% Node','Dimension','Value');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','dirichlet_data = [');
            fprintf(fileID,'\n');
            for k = 1:obj.nodes.vert
                fprintf(fileID,'%g %g %g\r\n',k,1,0);
                fprintf(fileID,'%g %g %g\r\n',k,2,0);
            end
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n');
            
            % Force presecribed
            fprintf(fileID,'%c','%% Force Prescribed');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s \t %s \t %s','% Node','Dimension','Value');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','pointload_complete = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n');
            
            % Volumetric force
            fprintf(fileID,'%c','%% Volumetric Force');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s \t %s \t %s','% Element','Dimension','Force_Dim');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','Vol_force = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n');
            
            % Group elements
            fprintf(fileID,'%c','%% Group Elements');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s \t %s','% Element','Group_num');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','Group = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n');
            
            % Initial holes
            fprintf(fileID,'%c','%% Group Elements');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','% Elements that are considered holes initially');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s','% Element');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','Initial_holes = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n');
            
            % Boundary elements
            fprintf(fileID,'%c','%% Boundary Elements');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','% Elements that can not be removed');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s','% Element');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','Boundary_elements = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n');
            
            % Boundary elements
            fprintf(fileID,'%c','%% Micro Gauss Post');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s','% Element');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','Micro_gauss_post = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n');
            
            % Master-slave nodes
            fprintf(fileID,'%c','%% Master-Slave');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','Master_slave = [');
            fprintf(fileID,'\n');
            for k = 1:size(obj.masterSlaveIndex,1)
                fprintf(fileID,'%g %g\r\n',obj.masterSlaveIndex(k,1),obj.masterSlaveIndex(k,2));
            end
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n');
            
            % Nodes solid
            fprintf(fileID,'%c','%% Nodes Solid');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','% Nodes that must remain');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s','% Nodes');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','% nodesolid = 1;');
            fprintf(fileID,'\n\n\n');
            
            % External border elements
            fprintf(fileID,'%c','%% External Border Elements');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','% Detect the elements that define the edge of the domain');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s \t %s \t %s','% Element','Node(1)','Node(2)');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','External_border_elements = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n');
            
            % External border nodes
            fprintf(fileID,'%c','%% External Border Nodes');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','% Detect the nodes that define the edge of the domain');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s \t %s \t %s','% Node');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','External_border_nodes = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
            fprintf(fileID,'\n\n\n');
            
            % Materials
            fprintf(fileID,'%c','%% Materials');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','% Materials that have been used');
            fprintf(fileID,'\n');
            fprintf(fileID,'%s \t %s \t %s \t %s','% Material_Num','Mat_density','Young_Modulus','Poisson');
            fprintf(fileID,'\n\n');
            fprintf(fileID,'%c','Materials = [');
            fprintf(fileID,'\n');
            fprintf(fileID,'%c','];');
        end
        
    end
    
end