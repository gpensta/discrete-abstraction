### discrete-abstraction

Add the CORA 2022 and functions directories to the path.

https://tumcps.github.io/CORA/pages/archive/index.html

Core files:

functions/get_alpha_max.m: Maximum spacing between target sets.

functions/get_ABW.m: First order linearization and computation of the linearization error set using Taylor models.

build_scc_graph.m: Constructs the graph of the discrete abstraction and colorizes the scc sets, Be aware that this file calls sccGraphviz.m, which in the last lines calls the Graphviz executable installed on Linux: adapt the command to your Graphviz installation.

demo_ref_control.m: Demonstrates how to determine if there is a path between two target sets.

nstep_perturbation.m: Computes the perturbation set when composing the inv pendulum dynamical functionmultiple times.