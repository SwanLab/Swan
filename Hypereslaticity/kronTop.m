function dom = kronTop(A,B)
    s.operation = @(xV) evaluate(A,B,xV);
    s.mesh      = A.mesh;
    s.ndimf     = A.ndimf;
    dom         = DomainFunction(s);
end

function C = evaluate(A,B,xV)
    Ae = A.evaluate(xV);
    Be = B.evaluate(xV);
    C = zeros([size(Ae,1), size(Ae,2), size(Be)]); % to support 4th order tensors
    for i = 1:size(Ae,1)
        for k = 1:size(Ae,2)
            for j = 1:size(Be,1)
                for l = 1:size(Be,2)
                    C(i,j,k,l,:,:) = Ae(i,k,:,:).*Be(j,l,:,:);
                end
            end
        end
    end
end