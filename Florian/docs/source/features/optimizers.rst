Alternative optimizers  
----------------------
    
The ``nullspace_optimizer`` package provides an interface with several optimizers   
which are available in the `nullspace_optimizer/optimizers <https://gitlab.com/florian.feppon/null-space-optimizer/-/tree/public-master/nullspace_optimizer/optimizers>`_   
folder. 
    
This interface allows to run these optimizers on the :class:`Optimizable` structure,
which is useful to compare the results with those of the function
:func:`~nullspace_optimizer.nlspace_solve`.     
    
Method of Moving Asymptotes 
^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
The function :func:`~nullspace_optimizer.optimizers.MMA.mma_solve` solves an
optimization problem using the implementation of the Method of Moving Asymptotes    
available on `this repository <https://github.com/arjendeetman/GCMMA-MMA-Python>`_.

.. code:: python
    
   from nullspace_optimizer.optimizers.MMA import mma_solve

.. autofunction:: nullspace_optimizer.optimizers.MMA.mma_solve
    
See `examples/topopt_examples/ex03_compliance_MMA.py <https://gitlab.com/florian.feppon/null-space-optimizer/-/blob/public-master/nullspace_optimizer//examples/topopt_examples/ex03_compliance_MMA.py>`_   
for an example of use.
    
.. note::   
    
   The decorator :func:`~nullspace_optimizer.bound_constraints_optimizable` should
   not be used when solving an optimization problem with  :func:`~nullspace_optimizer.optimizers.MMA.mma_solve`.
   Bound constraints need to be provided with the ``l`` and ``u`` arguments.

IPOPT   
^^^^^
    
The function :func:`~nullspace_optimizer.optimizers.IPOPT.ipopt_solve` solves an
optimization problem using the optimizer `IPOPT <https://coin-or.github.io/Ipopt/>`_ called by
the Python interface `cyipopt <https://cyipopt.readthedocs.io/en/stable/index.html>`_, which    
both need to be installed before running this function.

.. code:: python
    
   from nullspace_optimizer.optimizers.IPOPT import ipopt_solve

.. autofunction:: nullspace_optimizer.optimizers.IPOPT.ipopt_solve
    
See `examples/topopt_examples/ex04_compliance_IPOPT.py <https://gitlab.com/florian.feppon/null-space-optimizer/-/blob/public-master/nullspace_optimizer//examples/topopt_examples/ex04_compliance_IPOPT.py>`_   
for an example of use.

.. note::   
    
   The decorator :func:`~nullspace_optimizer.bound_constraints_optimizable` should
   not be used when solving an optimization problem with :func:`~nullspace_optimizer.optimizers.IPOPT.ipopt_solve`.     
   Bound constraints need to be provided with the ``l`` and ``u`` arguments.
    
Optimality Criteria 
^^^^^^^^^^^^^^^^^^^ 
    
The function :func:`~nullspace_optimizer.optimizers.OC.oc_solve` solves an
optimization problem using the Optimality Criteria method, based on its
implementation in the topology optimization code available `here <https://www.topopt.mek.dtu.dk/-/media/subsites/topopt/apps/dokumenter-og-filer-til-apps/topopt_cholmod.py>`_.
It can be used for solving optimization problems featuring a single equality
constraint and no inequality constraints apart from bound constraints. 

.. code:: python
    
   from nullspace_optimizer.optimizers.OC import oc_solve

.. autofunction:: nullspace_optimizer.optimizers.OC.oc_solve    
    
See `examples/topopt_examples/ex02_compliance_oc.py <https://gitlab.com/florian.feppon/null-space-optimizer/-/blob/public-master/nullspace_optimizer//examples/topopt_examples/ex02_compliance_oc.py>`_    
for an example of use.
    
.. note::   
    
   The decorator :func:`~nullspace_optimizer.bound_constraints_optimizable` should
   not be used when solving an optimization problem with  :func:`~nullspace_optimizer.optimizers.OC.oc_solve`.
   Bound constraints need to be provided with the ``l`` and ``u`` arguments.
