%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nonRegressionTest.m                           |
%|    #    |   VERSION    : 0.55                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : run all non regression test                   |
%|  `---'  |                                                              |
%+========================================================================+

clear all
close all
clc

% Gypsilab path
run('../addpathGypsilab.m')

% Mesh management
cd('meshManagement/')
run('nrtMshBitree.m')
run('nrtMshClean.m');
run('nrtMshCube.m');
run('nrtMshOctree.m')
run('nrtMshRead.m');
run('nrtMshRefine.m');
run('nrtMshSegment.m');
run('nrtMshSquare.m');
run('nrtMshTransfo.m');
run('nrtMshWrite.m');

% Domain quadrature
cd('../domainQuadrature/')
run('nrtDom1D.m');
run('nrtDomEdge.m');
run('nrtDomHalfSquare.m');
run('nrtDomIntSing2D.m');
run('nrtDomIntSing3D.m');
run('nrtDomND.m');
run('nrtDomTetra.m');
run('nrtDomTrace.m');
run('nrtDomTriangle.m');

% Finite element
cd('../finiteElement/')
run('nrtFemAllWithShuffle.m');
run('nrtFemContinuity1.m');
run('nrtFemContinuity2.m');
run('nrtFemConvergence2D.m');
run('nrtFemDirichletCube.m');
run('nrtFemDirichletDisk.m');
run('nrtFemDirichletSquare.m');
run('nrtFemDirichletString.m');
run('nrtFemElasticite2D')
run('nrtFemJunction.m');
run('nrtFemLaplace.m');
run('nrtFemOperators.m');
run('nrtFemRwgNed.m');
run('nrtFemWave1D.m');
run('nrtFemWave2D.m');

% Hierarchical matrices
cd('../hierarchicalMatrix')
run('nrtHmxAlgebra.m');
run('nrtHmxBEMConvergence.m');
run('nrtHmxBuilder.m');
run('nrtHmxBuilderFem.m');
run('nrtHmxCompareLU.m');
run('nrtHmxCompressorPartial.m');
run('nrtHmxCompressorsBox.m');
run('nrtHmxCompressorSingular');
run('nrtHmxCompressorTotal.m');
run('nrtHmxCriticalDimension.m');
run('nrtHmxCriticalForm.m');
run('nrtHmxLowrank.m');
run('nrtHmxNonSquare.m');

% Fast & Free memory Method
cd('../fastFreeMemory')
run('nrtFfmAlgebra.m');
run('nrtFfmBuilder.m');
run('nrtFfmBuilderFem.m');
run('nrtFfmHelmholtzBWdir.m');
run('nrtFfmHelmholtzBWneu.m');
run('nrtFfmMaxwellCFIEpec.m');

% Block matrix
cd('../blockMatrix')
run('nrtBmmAlgebra.m');
run('nrtBmmStokesConv.m');

% Scattering 2D
cd('../scattering2d')
run('nrtLaplace2dSDrad.m')
run('nrtHelmholtz2dSDrad.m')
run('nrtHmxHelmholtz2dS.m');
run('nrtHmxHelmholtz2dD.m');
run('nrtHmxHelmholtz2dDt.m');
run('nrtHmxHelmholtz2dH.m');
run('nrtHmxHelmholtz2dBWdir.m');
run('nrtHmxHelmholtz2dBWneu.m');

% Scattering 3D
cd('../scattering3d')
run('nrtHelmholtzCalderon.m')
run('nrtHelmholtzSDrad.m');
run('nrtHmxHelmholtzS.m');
run('nrtHmxHelmholtzD.m');
run('nrtHmxHelmholtzDxy.m');
run('nrtHmxHelmholtzDt.m');
run('nrtHmxHelmholtzH.m');
run('nrtHmxHelmholtzBWdir.m');
run('nrtHmxHelmholtzBWneu.m');
run('nrtHmxMaxwellT.m');
run('nrtHmxMaxwellNxK.m');
run('nrtHmxMaxwellCFIE.m');

% Fem-Bem dielectrique
cd('../femBemDielectrique')
run('nrtHmxFemBemEFIE.m');
run('nrtHmxFemBemEFIEhalf.m');
run('nrtHmxFemBemCFIE.m');
run('nrtHmxFemBemCFIEhalf.m');

% Inverse problem
cd('../inverseProblem')
run('nrtIpbHelmholtz.m');
run('nrtIpbHelmholtz0.m');

% Ray-tracing
cd('../rayTracing');
run('nrtRayCube.m');
run('nrtRayFabryPerot.m');
run('nrtRayLabyrinthe.m');
run('nrtRaySource.m');
run('nrtRaySphere.m');
run('nrtRayTheatre');

% Stokes
cd('../stokes');
run('nrtHmxStkConvergence.m');
run('nrtStkConvergence.m');
run('nrtStkRadiation.m');
run('translatingSphere/nrtStkTranslatingSphere.m');

% Vibro-acoustic
cd('../vibroAcoustic');
run('nrtVibroSlab2d.m');
run('nrtHmxVibroSlab2d.m');

% Operators
cd('../operators');
run('nrtOprValidation');

% End
disp('~~> Non Regresion Test are done. Michto gypsilab !');


