IMPORTANTE: modificada la clase lossFunctional mínimamente en las funciones computeStochasticCostAndGradient y computeFunctionAndGradient

De cara al debugging final en Julia: hacer ctrl+f ATTENTION para encontrar los transpuestos en la comunicacion intermadia (clases juliaCLASS.m)

Poner guion bajo delante de las propiedades y metodos privados.

SGD and Trainer have not been tested due to the added difficulty of having one inherit from the other. However they have been deeply examined to make them as easy to debug as possible once the whole code is run with Julia.

Write a test script for Nesterov, RMSprop and fminunc.

PlotterNN: uses getOutput() multiple times from costFunction but it's not defined. 
	Option 1) Define a getOutput() function in 		  costNN
	Option 2) See if it's just the name that is  		  outdated and then substitude 		  		  it by the name of the correct		  		  function.

Possible to change the remaining modules' constructors to the general nomenclature MODULEStruct