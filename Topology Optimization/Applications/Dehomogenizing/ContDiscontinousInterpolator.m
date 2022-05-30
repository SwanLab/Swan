classdef ContDiscontinousInterpolator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
       mesh
      % meshCont
      % fieldCont
       vector
       vectorCoherent

    end
    
    methods (Access = public)
        
        function obj = ContDiscontinousInterpolator(cParams)
            obj.init(cParams);            
        end
        
        function I = compute(obj)
            obj.computeCoherentVector()
            I = obj.computeI();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           % obj.meshCont = cParams.meshCont;
            obj.mesh  = cParams.meshDisc;
            obj.vector     = cParams.fieldDisc;
        end

        function computeCoherentVector(obj)
            nnode   = obj.mesh.nnodeElem;
            connec  = obj.mesh.connec;
            nElem = obj.mesh.nelem;            
            vD    = obj.vector;
            vC    = zeros(size(vD));
            for iElem = 1:nElem
                node1 = connec(iElem,1);
                v1 = vD(node1,:); 
                v2 = vD(node2,:); 
                v3 = vD(node3,:); 


                
                vC(node1,:) = v1;
                for iNode = 2:nnode
                    nodeI = connec(iElem,iNode);
                    vI = vD(nodeI,:);
                    z = cross([v1 0],[vI 0])/norm(v1)*norm(vI);
                    f1fI = -z(3);
                    vC(nodeI,:) = sign(f1fI)*vI;
                end                
            end
            obj.vectorCoherent = vC;
        end

        function I = computeI(obj)
            nnode   = obj.mesh.nnodeElem;
            connec = obj.mesh.connec;
            nElem = obj.mesh.nelem;
            vD    = obj.vector;
            vC    = obj.vectorCoherent;
            I = false(nnode,nnode,nElem);
            for iElem = 1:nElem
                for iNode = 1:nnode
                    nodeI = connec(iElem,iNode);
                    fI = vD(nodeI,:);
                    for jNode = 1:nnode
                        nodeJ = connec(iElem,jNode);
                        fJ = vC(nodeJ,:);
                        fIfJ = cos(dot(fI,fJ));
                        I(iNode,jNode,iElem) = sign(fIfJ);
                    end
                end
            end
        end
        
    end
    
end