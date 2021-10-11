%% FemTests
% Te sentit unificar, molts testos són identics i simplement en varia la
% geometria. Queda relativament petit, es poden agrupar tots els tipus de
% test en una mateixa classe.

%% UnfittedIntegrationTests
% Te prou sentit unificar, estructura pràcticament identica un cop s'agrupa
% tot. Nota: meshIncludeContour no fa absolutament res. Queda relativament
% petit, es poden agrupar tots els tipus de test en una mateixa classe.

%% VectorizedTriangulationTests o CutCellsTests
% Aqui hi ha molt de joc amb el 'Sequential', pero el Bcutmesh pot ser
% bastant plasta. Una opcio es fer-los tots amb i sense, tampoc hauria de
% venir d'aqui. El mateix amb els randoms.

%% TopOptTests
% Bona part son agrupables, n'hi ha que no, n'hi ha que no funcionen

%% ReadingFilesTests
% Nomes n'hi ha un :( Pero es pot reciclar de FemTests el
% PrecomputedVariableTest

%% PlottingTests
% S'hi pot fer feina 

%% HomogenizationTests
% Són molt diversos i es fa dificil crear una arquitectura decent per
% englobar-ne a més d un. El que tindria més sentit és canviar els
% testShowingError per UnitTests però matenir-ho tot separat. Que la Suite
% faci correr tots els testos, en classes diferents

%% ImageProcessingTests
% Es pot fer, i es pot reciclar parcialment codi de FemTests

FemTests;
UnfittedIntegrationTests;
VectorizedTriangulationTests;
TopOptTests; % 2 tests amb problemes
ReadingFilesTests;
PlottingTests;
HomogenizationTests;
ImageProcessingTests;