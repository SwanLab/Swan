In this folder, you can find five main functions:

    - testMULTIPROD.m
    - timing_MX.m
    - timing_matlab_commands.m
    - timing_arraylab_engines.m
    - testing_memory_usage.m

testMULTIPROD checks the results of 15979 different kinds of multiproducts,
most of which require array expansion (see also Appendix C).

The other four main functions were used to test the speed and memory usage
of MULTIPROD, its engines and some MATLAB commands. This allowed us to take
decisions regarding the optimization of MULTIPROD (see Appendix B).

This folder also contains a few secondary functions, which are just meant
to be called by the above-listed main functions.