function dom = kronProd(A,B,idx)
    s.operation = @(xV) evaluate(A,B,idx,xV);
    s.mesh      = A.mesh;
    s.ndimf     = A.ndimf;
    dom         = DomainFunction(s);
end

function C = evaluate(A,B,idx,xV)
    Ae = Expand(A,2).evaluate(xV);
    Be = Expand(B,2).evaluate(xV);
    C = zeros([size(Ae,[1 2]), size(Be)]);
    for i = 1:size(Ae,1)
        for j = 1:size(Ae,2)
            for k = 1:size(Be,1)
                for l = 1:size(Be,2)
                    C(i,j,k,l,:,:) = Ae(i,j,:,:).*Be(k,l,:,:);
                end
            end
        end
    end
    C = squeezeParticular(C,[1 2 3 4]);
    nF = length(idx);
    [~,pos] = sort(idx);
    C = permute(C,[pos nF+1 nF+2]);
end