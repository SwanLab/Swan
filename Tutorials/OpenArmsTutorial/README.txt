IMPORTANTE: modificada la clase lossFunctional mínimamente en las funciones computeStochasticCostAndGradient y computeFunctionAndGradient

De cara al debugging final en Julia: hacer ctrl+f ATTENTION para encontrar los transpuestos en la comunicacion intermadia (clases juliaCLASS.m)

Poner guion bajo delante de las propiedades y metodos privados.

SGD and Trainer have not been tested due to the added difficulty of having one inherit from the other. However they have been deeply examined to make them as easy to debug as possible once the whole code is run with Julia.