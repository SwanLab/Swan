function fVR = Integral(operation,mesh,quad)
    xV = quad.posgp;
    dVolu(1,1,:,:)  = mesh.computeDvolume(quad);
    integ = bsxfun(@times, operation.evaluate(xV), dVolu);
    integ = squeezeParticular(sum(integ,3),3);
    fVR = integ;
end

% function dom = Integral(operation,mesh)
%     s.operation = @(quad) evaluate(operation, mesh, quad);
%     dom = IntegralFunction(s);
% end
% 
% function fVR = evaluate(operation, mesh, quad)
%     xV = quad.posgp;
%     dVolu(1,1,:,:)  = mesh.computeDvolume(quad);
%     integ = bsxfun(@times, operation.evaluate(xV), dVolu);
%     integ = squeezeParticular(sum(integ,3),3);
%     fVR = integ;
% end