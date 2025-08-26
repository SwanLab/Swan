function dom = kronProd(A,B,idx)
    s.operation = @(xV) evaluate(A,B,idx,xV);
    s.mesh      = A.mesh;
    s.ndimf     = A.ndimf;
    dom         = DomainFunction(s);
end

function C = evaluate(A,B,idx,xV)
    Ae = A.evaluate(xV);
    Be = B.evaluate(xV);
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
    [~,pos] = sort(idx);
    C = permute(C,[pos 5 6]);
end