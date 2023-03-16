% eg. RaviartThomasElement2D -> create public computeShapeDeriv(obj) method
%   crear propietat publica "shape" amb les funcions de forma
%   no passa res per no tenir les derivades, amb les shape n'hi ha prou

% crear-se FeLagrangianFunction tipus P1function mes o menys
% la clau estara en crear-se be la interpolation. fins ara, tenim
% d'interpolation coses com Triangle_Linear, s'ha de substituir pels teus
% elements

% dins de la FeLagrangianFunction, crear-se el metode evaluate, que fara
% servir la teva interpolation