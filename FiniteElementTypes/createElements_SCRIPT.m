%% LAGRANGE SIMPLICIAL
k = 1;
d = 1;
ls = LagrangeElement.create('SIMPLICIAL',k,d);
ls.plotshapefunctions();

%% LAGRANGE TENSOR PRODUCT
k = 1;
d = 1;
lt = LagrangeElement.create('TENSOR PRODUCT',k,d);
lt.plotshapefunctions();

%% CROUZEIX-RAVIART
d = 2;
cr = CrouzeixRaviartElement.create(d);
cr.plotshapefunctions();

%% RAVIART-THOMAS
d = 2;
rt = RaviartThomasElement.create(d);
rtplotshapefunctions();

%% NÉDÉLEC
d = 2;
nd = NedelecElement.create(d);
nd.plotshapefunctions();
