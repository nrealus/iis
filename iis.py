"""
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

Since the optimization problems that appear in [2] are NLPs (where only one
constraint is non linear), the MIP obtained by piecewise linearization isn't
a "real" MIP in the sense that the integer / binary variables are auxiliary
variables used for the discretization / linearization (SOS2 - type 2 Special
Ordered Sets variables), and supported internally by most LP/MIP solvers
(including CBC, the default one interfaced with Python-MIP). So it actually
wouldn't be informative to have these auxiliary linearization constraints as
part of the IISs. (FIXME: Still unsure about that though...) As such, the
fact that the SOS2 constraints are "invisible" for us from the Python-MIP API
(and therefore are not included in the IIS we return) should not be seen as a problem, in our case.
(FIXME: again, unsure).

[1]: John W. Chinneck. Feasibility and Infeasibility in Optimization: Algorithms and Computational Methods (2008)
[2]: Wang, A.J. Risk-Bounded Dynamic Scheduling of Temporal Plans (2022)

"""

# TODO more tests

from __future__ import annotations

from typing import Optional, Set

import mip

#################################################################################


#################################################################################
#
#################################################################################
 
def get_iis_additive_deletion_method(
    mip_model: mip.Model,
    premade_aux_mip_model: Optional[mip.Model] = None,
) -> mip.ConstrList:
    
    # TODO: Add ability to specify constraints that should be
    # ignored / excluded from the returned IIS ? The motivation is
    # that we may not want to have "auxiliary" constraints in the IIS.
    # Indeed, intuitively, some constraints (for example those introduced
    # for piecewise linearization using SOS2) are semantically "bundled",
    # as they participate to express one "macro" constraint.
    #
    # (does this even make sense ?... lol)

    if mip_model.status == mip.OptimizationStatus.LOADED:
        mip_model.optimize()

    if (mip_model.status == mip.OptimizationStatus.FEASIBLE
        or mip_model.status == mip.OptimizationStatus.OPTIMAL
    ):
        return set()    # type: ignore

    if premade_aux_mip_model is not None:
        aux_mip_model = premade_aux_mip_model
        aux_mip_model.clear()
    else:
        aux_mip_model = mip.Model()
        aux_mip_model.verbose = 0

    aux_mip_model.emphasis = mip.SearchEmphasis.FEASIBILITY
    aux_mip_model.preprocess = 1

    aux_mip_model.objective = 0     # type: ignore

    for var in mip_model.vars:
        aux_mip_model.add_var(name=var.name,
                              lb=var.lb,
                              ub=var.ub,
                              var_type=var.var_type)

    iis = aux_mip_model.constrs
    
    i = 0
    for constr in mip_model.constrs:

        iis.add(constr.expr)
        aux_mip_model.optimize()

        if (aux_mip_model.status == mip.OptimizationStatus.INFEASIBLE
            or aux_mip_model.status == mip.OptimizationStatus.INT_INFEASIBLE
        ):
            break
        i += 1

    temp = iis[i:].copy()       # type: ignore
#    for constr in iis[i:]:      # type: ignore
    for constr in temp:         # type: ignore
        expr = constr.expr

        iis.remove([constr])
        aux_mip_model.optimize()

        if (aux_mip_model.status == mip.OptimizationStatus.FEASIBLE
            or aux_mip_model.status == mip.OptimizationStatus.OPTIMAL
        ):
            aux_mip_model.add_constr(expr)     
            # /!\ using constr.expr instead of "cached" expr results
            # in an "invalid row index (-1) ..." error from CBC. This is
            # probably due to internal logic in "ConstrList.remove"
            # (see "iis.remove([constr])" above).
    
    return iis

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest
import statistics

class Tests(unittest.TestCase):

    def assertConstrExprSetEqual(self,
        constr_exprs_set1: Set[mip.LinExpr],
        constr_exprs_set2: Set[mip.LinExpr],
    ):
        ok = True

        for lin_expr1 in constr_exprs_set1:
            ok = False
            for lin_expr2 in constr_exprs_set2:

                if lin_expr1.equals(lin_expr2):
                    ok = True
                    break

            if ok == False:
                msg = "Sets of mip.LinExpr are not identical."
                raise self.failureException(msg)

        for lin_expr2 in constr_exprs_set2:
            ok = False
            for lin_expr1 in constr_exprs_set1:

                if lin_expr1.equals(lin_expr2):
                    ok = True
                    break

            if ok == False:
                msg = "Sets of mip.LinExpr are not identical."
                raise self.failureException(msg)

    def test01_iis(self):

        m = mip.Model()
        m.verbose = 0

        v1 = m.add_var()
        v2 = m.add_var()

        c1 = m.add_constr(v1 <= 5)      # type: ignore
        c2 = m.add_constr(v1 >= 6)      # type: ignore
        m.add_constr(v2 <= 1)           # type: ignore

        m.optimize()

        self.assertTrue(m.status == mip.OptimizationStatus.INFEASIBLE)

        iis = get_iis_additive_deletion_method(m)
 
        self.assertConstrExprSetEqual({c.expr for c in iis},
                                      set((c1.expr, c2.expr)))

    def test02_iis(self):

        mu = 4
        sigma = 1
        a = 1
        b = 7

        num_discr_pts = 10
        step = (b - a) / num_discr_pts
        discr_pts = tuple(a + i*step for i in range(num_discr_pts))

        normal = statistics.NormalDist(0, 1)
        alpha = (a - mu) / sigma
        beta = (b - mu) / sigma

        def truncated_normal_cdf(x):
            return (normal.cdf((x - mu) / sigma) - normal.cdf(alpha)) / (normal.cdf(beta) - normal.cdf(alpha))
        
        m = mip.Model()
        m.verbose = 0

        l = m.add_var()
        F_l = m.add_var()

        u = m.add_var()
        F_u = m.add_var()

        w_l = [m.add_var(lb=0, ub=1, var_type=mip.CONTINUOUS) for _ in range(num_discr_pts)]    # type: ignore
        m.add_sos([(w_l[k], discr_pts[k]) for k in range(num_discr_pts)], 2)                    # type: ignore
        m.add_constr(mip.xsum(w_l) == 1)
        m.add_constr(l == mip.xsum(discr_pts[k] * w_l[k] for k in range(num_discr_pts)))        # type: ignore
        m.add_constr(F_l == mip.xsum(truncated_normal_cdf(discr_pts[k]) * w_l[k] for k in range(num_discr_pts)))        # type: ignore

        w_u = [m.add_var(lb=0, ub=1, var_type=mip.CONTINUOUS) for _ in range(num_discr_pts)]    # type: ignore
        m.add_sos([(w_u[k], discr_pts[k]) for k in range(num_discr_pts)], 2)                    # type: ignore
        m.add_constr(mip.xsum(w_u) == 1)
        m.add_constr(u == mip.xsum(discr_pts[k] * w_u[k] for k in range(num_discr_pts)))        # type: ignore
        m.add_constr(F_u == mip.xsum(truncated_normal_cdf(discr_pts[k]) * w_u[k] for k in range(num_discr_pts)))        # type: ignore

        m.add_constr(F_u - F_l >= 0.8)  # type: ignore
        m.add_constr(u - l <= 2)        # type: ignore

        m.optimize()

        self.assertTrue(m.status == mip.OptimizationStatus.INFEASIBLE)

        iis = get_iis_additive_deletion_method(m)

        self.assertConstrExprSetEqual({c.expr for c in iis},
                                      set((c.expr for c in m.constrs)))

    def test03_iis(self):

        mu = 4
        sigma = 1
        a = 1
        b = 7

        num_discr_pts = 10
        step = (b - a) / num_discr_pts
        discr_pts = tuple(a + i*step for i in range(num_discr_pts))

        normal = statistics.NormalDist(0, 1)
        alpha = (a - mu) / sigma
        beta = (b - mu) / sigma

        def truncated_normal_cdf(x):
            return (normal.cdf((x - mu) / sigma) - normal.cdf(alpha)) / (normal.cdf(beta) - normal.cdf(alpha))
        
        m = mip.Model()
        m.verbose = 0

#        l = m.add_var()
#        F_l = m.add_var()

#        u = m.add_var()
#        F_u = m.add_var()

        w_l = [m.add_var(lb=0, ub=1, var_type=mip.CONTINUOUS) for _ in range(num_discr_pts)]    # type: ignore
        m.add_sos([(w_l[k], discr_pts[k]) for k in range(num_discr_pts)], 2)                    # type: ignore
        m.add_constr(mip.xsum(w_l) == 1)
#        m.add_constr(l == mip.xsum(discr_pts[k] * w_l[k] for k in range(num_discr_pts)))           # type: ignore
        l = mip.xsum(discr_pts[k] * w_l[k] for k in range(num_discr_pts))                           # type: ignore
#        m.add_constr(F_l == mip.xsum(truncated_normal_cdf(discr_pts[k]) * w_l[k] for k in range(num_discr_pts)))        # type: ignore
        F_l = mip.xsum(truncated_normal_cdf(discr_pts[k]) * w_l[k] for k in range(num_discr_pts))   # type: ignore

        w_u = [m.add_var(lb=0, ub=1, var_type=mip.CONTINUOUS) for _ in range(num_discr_pts)]    # type: ignore
        m.add_sos([(w_u[k], discr_pts[k]) for k in range(num_discr_pts)], 2)                    # type: ignore
        m.add_constr(mip.xsum(w_u) == 1)
#        m.add_constr(u == mip.xsum(discr_pts[k] * w_u[k] for k in range(num_discr_pts)))       # type: ignore
        u = mip.xsum(discr_pts[k] * w_u[k] for k in range(num_discr_pts))                       # type: ignore
#        m.add_constr(F_u == mip.xsum(truncated_normal_cdf(discr_pts[k]) * w_u[k] for k in range(num_discr_pts)))        # type: ignore
        F_u = mip.xsum(truncated_normal_cdf(discr_pts[k]) * w_u[k] for k in range(num_discr_pts))   # type: ignore

        m.add_constr(F_u - F_l <= 0.3)  # type: ignore
        m.add_constr(u - l >= 4)  # type: ignore

        m.optimize()

        self.assertTrue(m.status == mip.OptimizationStatus.INFEASIBLE)

        iis = get_iis_additive_deletion_method(m)

        self.assertConstrExprSetEqual({c.expr for c in iis},
                                      set((c.expr for c in m.constrs)))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
    unittest.main()