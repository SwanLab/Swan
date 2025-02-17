%MAIN for TestNaca

Naca.M     = 0.02;
Naca.p     = 0.4;
Naca.t     = 12/100;
Naca.chord = 1;
Naca.AoA   = 0;

NacaClass = TestNaca(Naca);
NacaClass.compute();
NacaClass.validate();
NacaClass.print();
