# Copyright 2018-2019 CNRS, Ecole Polytechnique and Safran.
#
# This file is part of nullspace_optimizer.
#
# nullspace_optimizer is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# nullspace_optimizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is included below.
# For further information, see <http://www.gnu.org/licenses/>.

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nullspace_optimizer",
    version="1.2.6",
    author="Florian Feppon, Dries Toebat",
    author_email="florian.feppon@kuleuven.be",
    license="GNU GPL version 3",
    description="Null space algorithm for nonlinear constrained optimization",
    keywords="nonlinear constrained optimization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://null-space-optimizer.readthedocs.io/en/latest/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Mathematics"
    ],
    install_requires=['numpy>=1.26.4',
                      'scipy>=1.11.4',
                      'matplotlib>=3.8.3',
                      'cvxopt>=1.3.2',
                      'osqp>=0.6.2',    
                      'sympy>=1.12',
                      'qpalm>=1.2.2',   
                      'colored>=1.4.4'
                      ],
    python_requires='>=3.9'
)
