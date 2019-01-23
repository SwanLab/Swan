function [kappa_opt] = kappa_imposed(element,iter)

hkappa = element.material.hkappa;
kappa_opt = hkappa(1,iter);

end
