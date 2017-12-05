classdef Postprocess
    properties
        Values
    end
    
    methods (Access = protected)
        function [nnode,ndim,pdim,gtype,etype,nelem,npnod,coordinates,conectivities] = getBasicParams(obj,input)
            coordinates = input.mesh.coord;
            conectivities = input.mesh.connec;
            gtype = input.mesh.geometryType;
            nnode = length(conectivities(1,:));
            ndim = input.dim.ndim;
            pdim = input.mesh.pdim;
            
            switch  gtype %gid type
                case 'TRIANGLE'
                    etype = 'Triangle';
                case 'QUAD'
                    etype = 'Quadrilateral';
                case 'TETRAHEDRA'
                    etype = 'Tetrahedra';
                case 'HEXAHEDRA'
                    etype = 'Hexahedra';
            end
            nelem = size(conectivities,1); % Number of elements
            npnod = size(coordinates,1);   % Number of nodes
        end
        function printTitle(obj,fid)
            fprintf(fid,'####################################################\n');
            fprintf(fid,'################# FEM-MAT-OO v.1.0 #################\n');
            fprintf(fid,'####################################################\n');
            fprintf(fid,'\n');
        end
        function PrintVector(obj,fid,ndim,nameres,indexName,problemType,result_type,istep,result_location,location_name,results)
            % Print Header ------------------------------------------------
            fprintf(fid,'\nResult "%s" "%s" %.0f %s %s "%s"\n',nameres,problemType,istep,result_type,result_location,location_name);
            switch ndim
                case 2
                    fprintf(fid,'ComponentNames  "%sx", "%sy"\n',indexName,indexName);
                case 3
                    fprintf(fid,'ComponentNames "%sx", "%sy", "%sz"\n',indexName,indexName,indexName);
                otherwise
                    error('Invalid value of parametre ndime.')
            end
            
            % Print Variables ---------------------------------------------
            fprintf(fid,'Values\n');
            for inode = 1:round(length(results)/ndim)
                fprintf(fid,'%6.0f ',inode);
                for idime = 1:ndim
                    fprintf(fid,'%12.5d ',results(ndim*(inode-1)+idime));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'End Values\n');
        end
        function PrintTensor(obj,fid,ndim,nameres,indexName,problemType,result_type,istep,result_location,location_name,results)
            % Print Header ------------------------------------------------
            fprintf(fid,'\nResult "%s" "%s" %.0f %s %s "%s"\n',nameres,problemType,istep,result_type,result_location,location_name);
            switch ndim
                case 2
                    fprintf(fid,'ComponentNames  "%sx", "%sy", "%sxy", "%sz"\n',indexName,indexName,indexName,indexName);
                case 3
                    fprintf(fid,'ComponentNames "%sx", "%sy", "%sz", "%sxy", "%syz", "%sxz"\n',indexName,indexName,indexName,indexName);
                otherwise
                    error('Invalid value of parametre ndime.')
            end
            
            % Print Variables ---------------------------------------------
            fprintf(fid,'Values\n');
            for ielem = 1:size(results,3)
                fprintf(fid,'%6.0f ',ielem);
                for igaus = 1:size(results,1)
                    for istre = 1:size(results,2)
                        fprintf(fid,'%12.5d ',results(igaus,istre,ielem));
                    end
                    fprintf(fid,'\n');
                end
            end
            fprintf(fid,'End Values\n');
        end
        function PrintEscalar(obj,fid,ndim,nameres,indexName,problemType,result_type,istep,result_location,location_name,results)
            % Print Header ------------------------------------------------
            fprintf(fid,'\nResult "%s" "%s" %.0f %s %s "%s"\n',nameres,problemType,istep,result_type,result_location,location_name);
            fprintf(fid,'ComponentNames  "%s"\n',indexName);                      
            
            % Print Variables ---------------------------------------------
            fprintf(fid,'Values\n');
            for inode = 1:length(results)
                fprintf(fid,'%6.0f ',inode);
                fprintf(fid,'%12.5d ',results(inode));
                fprintf(fid,'\n');
            end            
            fprintf(fid,'End Values\n');
        end
    end
end

