# Additive-deletion filter method for IIS computation with Python-MIP

The function implemented in this module can be used to get the IIS (Irreducible
Inconsistent Set (or Subsystem)) of an infeasible LP or MIP, in the Python-MIP library.
But the method used can also be applied to NLPs.

The method used is one of the simpler algorithms for IIS computation, called
the "additive/deletion filter" [1] (Chapter 6.1.6 / Algorithm 6.7).

Many of the algorithms (including the method used) presented in that book have
specialized versions for the cases of LPs or MIPs, for more efficiency.
However, here we only implement the basic general method.

This code was developed as part of an effort to implement the algorithm for
risk-bounded dynamic controllabilty checking of PSTNs [2], with some optimizations,
namely piecewise linearization of NLP, and use of IIS as the set of conflicting
constraints instead of all of them.

[1]: John W. Chinneck. Feasibility and Infeasibility in Optimization: Algorithms and Computational Methods (2008)
[2]: Wang, A.J. Risk-Bounded Dynamic Scheduling of Temporal Plans (2022)
