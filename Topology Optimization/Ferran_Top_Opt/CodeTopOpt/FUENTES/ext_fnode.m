function [extfnod] = ext_fnode(nndof,nunkn,pointload,fvol,ftype)

switch ftype
    case {'ELASTIC'}
        nodes=size(pointload,1);
        extfnod = zeros(nndof,1);
        for i=1:nodes
            ipoin = pointload(i,1);
            idime = pointload(i,2);
            iglib = nunkn*(ipoin-1)+idime;
            extfnod(iglib)=extfnod(iglib)+pointload(i,3);
        end
        
        
    case {'THERMAL'}
        nodes=size(pointload,1);
        extfnod = zeros(nndof,1);
        %in thermal x direction refers to normal direction of the load.
        for i=1:nodes
            if pointload(i,2) == 1
            ipoin = pointload(i,1);
            idime = pointload(i,2);
            iglib = nunkn*(ipoin-1)+idime;
            extfnod(iglib)=extfnod(iglib)+pointload(i,3);
            end
        end

        
end

extfnod = extfnod + fvol;
end

