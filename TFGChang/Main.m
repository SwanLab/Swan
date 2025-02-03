%MAIN for TestNaca

Naca.p = 0.5;
Naca.M = 0.02;
Naca.t = 12/100;

NacaClass = TestNaca(Naca);
NacaClass.compute();
NacaClass.validate();