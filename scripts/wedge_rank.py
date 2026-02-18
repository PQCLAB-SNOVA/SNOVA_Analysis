# https://gitlab.di.ens.fr/hle/snova_wedge_attack
# commit 87f695020642701e26c5109666a9c8d06fcecf83
"""
Rank predictions for SNOVA wedge attack.

We have two models:
  - Using m 2-forms | Scenario 1 | Proposition 2 | Conservative:  N^{<1>}  => cross term exponent 1  ( (1 - t*y*z_i z_j)^(-1) )^m
  - Using ml 2-forms | Scenario 2 | Proposition 3 | Aggressive:    N^{<2>}  => cross term exponent 2  ( (1 - t*y*z_i z_j)^(-2) )^m

Rankpred is (no min cap of number of variables):
  rankpred = S

This script supports:
  - General (n, v, m, l)
  - Project ANY subset of the l blocks and choose kept oil per block
  - Compute (S, vars) for conservative/aggressive
  - Search for best projection under two projection regimes:
      * unbalanced projection (default): general sorted tuples (multiset search)
      * balanced projection: tuples with equal elements (x, x, ..., x)

Block model:
  - There are l blocks; each block i keeps o'_i oil coords (0..o=n-v),
    thus block size n_i = v + o'_i.
  - vars = Π_i binom(n_i, v)

Usage examples are at the bottom.
"""

from __future__ import annotations
from collections import defaultdict
from functools import lru_cache
from math import comb
import math
import heapq
from typing import Dict, Iterable, List, Optional, Tuple, Union

# ---------------------------------------------------------------------
# Projection helpers
# ---------------------------------------------------------------------


def comb0(n: int, k: int) -> int:
    """math.comb but returns 0 when k is out of range."""
    return comb(n, k) if 0 <= k <= n else 0


def compute_n_eqs(v: int, m: int, l: int, n_blocks: List[int], mode: str) -> Tuple[int, int]:
    """
    n_eqs = sum_{a_0+...+a_{l-1}=2} M * Π_j C(n_blocks[j], v + a_j)

    Interprets the "v+o'" per block as n_blocks[j] = v + o'_j.
    M = m (conservative), M = m*l (aggressive).
    """
    if mode not in ("conservative", "aggressive"):
        raise ValueError("mode must be 'conservative' or 'aggressive'")

    M = m if mode == "conservative" else m * l

    B0 = [comb0(ni, v) for ni in n_blocks]       # a_j = 0
    B1 = [comb0(ni, v + 1) for ni in n_blocks]   # a_j = 1
    B2 = [comb0(ni, v + 2) for ni in n_blocks]   # a_j = 2

    P0 = 1
    for x in B0:
        P0 *= x  # Π_j C(n_j, v)

    # sum over patterns with total degree 2:
    # (i) one block gets 2, others 0
    term_one_2 = 0
    for j in range(l):
        if B2[j] == 0:
            continue
        term_one_2 += B2[j] * (P0 // B0[j])

    # (ii) two blocks get 1,1; others 0
    term_two_1 = 0
    for i in range(l):
        if B1[i] == 0:
            continue
        for j in range(i + 1, l):
            if B1[j] == 0:
                continue
            denom = B0[i] * B0[j]
            term_two_1 += B1[i] * B1[j] * (P0 // denom)

    n_eqs = M * (term_one_2 + term_two_1)
    return n_eqs, M


def make_o_kept_vector(
    n: int,
    v: int,
    l: int,
    proj_blocks: Optional[Iterable[int]] = None,
    o_kept: Union[int, List[int], Tuple[int, ...], Dict[int, int], None] = None,
) -> List[int]:
    """
    Return o_kept_vec of length l with 0 <= o'_i <= o=n-v.

    - If proj_blocks is None: no projection => all blocks keep full o.
    - Else proj_blocks specifies which blocks are projected.
      Blocks not in proj_blocks keep full o.
      o_kept can be:
        * int => same o' for all projected blocks
        * list/tuple (len = len(proj_blocks)) => per projected block
        * dict {block_index: o'} => explicit per projected block
    """
    o = n - v
    if proj_blocks is None:
        return [o] * l

    proj_blocks = list(proj_blocks)
    if any((b < 0 or b >= l) for b in proj_blocks):
        raise ValueError(f"proj_blocks must be in [0..{l-1}]")

    out = [o] * l  # default: keep full oil

    if isinstance(o_kept, int):
        for b in proj_blocks:
            out[b] = o_kept
    elif isinstance(o_kept, dict):
        for b in proj_blocks:
            if b not in o_kept:
                raise ValueError(f"Missing o_kept for projected block {b} in dict.")
            out[b] = o_kept[b]
    elif isinstance(o_kept, (list, tuple)):
        if len(o_kept) != len(proj_blocks):
            raise ValueError("If o_kept is list/tuple, it must match len(proj_blocks).")
        for b, oi in zip(proj_blocks, o_kept):
            out[b] = oi
    else:
        raise ValueError("o_kept must be int, dict, list/tuple, or None.")

    for i, oi in enumerate(out):
        if not (0 <= oi <= o):
            raise ValueError(f"Block {i}: o_kept={oi} out of range [0, {o}].")
    return out


# ---------------------------------------------------------------------
# Generating-series coefficient engine
# ---------------------------------------------------------------------

def _convolve_state(
    state: Dict[Tuple[int, Tuple[int, ...]], int],
    terms: List[Tuple[int, Dict[int, int], int]],
    o_bounds: List[int],
    rmax: int,
) -> Dict[Tuple[int, Tuple[int, ...]], int]:
    """
    state[(r, a_tuple)] = coeff
    terms entries are (dr, {idx: delta_ai}, coeff)
    """
    new_state = defaultdict(int)
    for (r0, a0), c0 in state.items():
        for dr, da, c1 in terms:
            r = r0 + dr
            if r > rmax:
                continue
            a = list(a0)
            ok = True
            for idx, inc in da.items():
                a[idx] += inc
                if a[idx] > o_bounds[idx]:
                    ok = False
                    break
            if not ok:
                continue
            new_state[(r, tuple(a))] += c0 * c1
    return dict(new_state)


@lru_cache(maxsize=None)
def _coeffs_Nr_sumk_cached(
    m: int,
    l: int,
    o_bounds_tup: Tuple[int, ...],
    cross_mult: int,
) -> Tuple[Dict[Tuple[int, Tuple[int, ...]], int], int]:
    """
    Cached wrapper so repeated search calls are cheaper.
    Returns (state, rmax) where state[(r,a)] = coefficient N_r(a).
    """
    o_bounds = list(o_bounds_tup)
    rmax = sum(o_bounds) // 2
    zero = tuple([0] * l)
    state: Dict[Tuple[int, Tuple[int, ...]], int] = {(0, zero): 1}

    # Self factors: (1 - t z_i^2)^(-m)
    # = Σ_s binom(m+s-1, s) t^s z_i^(2s)
    for i in range(l):
        max_s = min(rmax, o_bounds[i] // 2)
        terms = []
        for s in range(max_s + 1):
            terms.append((s, {i: 2 * s}, comb(m + s - 1, s)))
        state = _convolve_state(state, terms, o_bounds, rmax)

    # Cross factors: (1 - t z_i z_j)^(-cross_mult*m)
    # = Σ_s binom(cross_mult*m + s - 1, s) t^s z_i^s z_j^s
    exp = cross_mult * m
    for i in range(l):
        for j in range(i + 1, l):
            max_s = min(rmax, o_bounds[i], o_bounds[j])
            terms = []
            for s in range(max_s + 1):
                terms.append((s, {i: s, j: s}, comb(exp + s - 1, s)))
            state = _convolve_state(state, terms, o_bounds, rmax)

    return state, rmax


def coeffs_Nr_sumk(m: int, l: int, o_bounds: List[int], cross_mult: int):
    """
    cross_mult = 1 => conservative N^{<1>}
    cross_mult = 2 => aggressive   N^{<2>}
    """
    return _coeffs_Nr_sumk_cached(m, l, tuple(o_bounds), cross_mult)


# ---------------------------------------------------------------------
# Compute S and vars for conservative/aggressive with per-block projection
# ---------------------------------------------------------------------

def compute_S_and_vars(
    n: int,
    v: int,
    m: int,
    l: int,
    o_kept_vec: List[int],
    mode: str,
) -> Dict[str, object]:
    """
    mode:
      - "conservative" => cross_mult = 1
      - "aggressive"   => cross_mult = 2

    rankpred = S (per Proposition 2 and 3)
    """
    if mode not in ("conservative", "aggressive"):
        raise ValueError("mode must be 'conservative' or 'aggressive'")

    o = n - v
    if len(o_kept_vec) != l:
        raise ValueError("o_kept_vec length must equal l")
    for oi in o_kept_vec:
        if not (0 <= oi <= o):
            raise ValueError(f"each o'_i must be in [0, {o}]")

    cross_mult = 1 if mode == "conservative" else 2

    n_blocks = [v + oi for oi in o_kept_vec]   # n_i
    o_bounds = list(o_kept_vec)                # bounds for a_i
    state, rmax = coeffs_Nr_sumk(m, l, o_bounds, cross_mult=cross_mult)

    # vars = Π_i binom(n_i, v)
    vars_ = 1
    for ni in n_blocks:
        vars_ *= comb(ni, v)

    # Number of eqs generated: n_eqs
    n_eqs, M = compute_n_eqs(v=v, m=m, l=l, n_blocks=n_blocks, mode=mode)

    # S = alternating sum over coefficients
    S = 0
    for (r, a), Nr in state.items():
        if r == 0 or Nr == 0:
            continue
        if sum(a) != 2 * r:
            continue
        prod = 1
        for i in range(l):
            prod *= comb(n_blocks[i], v + a[i])
        sign = 1 if (r % 2 == 1) else -1  # (-1)^(r+1)
        S += sign * Nr * prod

    return {
        "mode": mode,
        "S": S,
        "vars": vars_,
        "rankpred": S,
        "n_eqs": n_eqs,
        "M": M,
        "o_kept_vec": list(o_kept_vec),
        "n_blocks": n_blocks,
        "rmax": rmax,
        "cross_mult": cross_mult,
    }


def compute_S_and_vars_from_projection(
    n: int,
    v: int,
    m: int,
    l: int,
    mode: str,
    proj_blocks: Optional[Iterable[int]] = None,
    o_kept: Union[int, List[int], Tuple[int, ...], Dict[int, int], None] = None,
) -> Dict[str, object]:
    """
    Convenience wrapper: specify projection as (proj_blocks, o_kept).
    """
    o_kept_vec = make_o_kept_vector(n, v, l, proj_blocks, o_kept)
    return compute_S_and_vars(n, v, m, l, o_kept_vec, mode=mode)


# ---------------------------------------------------------------------
# Search: minimize vars subject to rankpred=S exceeding vars
# ---------------------------------------------------------------------

def find_best_projection(
    n: int,
    v: int,
    m: int,
    l: int,
    mode: str,
    projection: str = "unbalanced",
    strict: bool = True,
    max_states: int = 20000,
) -> Dict[str, object]:
    """
    Minimize vars(o') subject to S(o') > vars(o') (strict) or >= (non-strict),
    where S is computed for the chosen mode.

    projection:
      - "unbalanced" (default): searches only sorted o'_vectors (multisets)
                                using best-first expansion by vars.
      - "balanced": restricts to tuples with equal elements, i.e. o'_1=...=o'_l,
                    scanning in increasing vars (equivalently increasing o').

    Returns the first feasible solution found (optimal for vars under the
    respective projection constraint).
    """
    if mode not in ("conservative", "aggressive"):
        raise ValueError("mode must be 'conservative' or 'aggressive'")
    if projection not in ("unbalanced", "balanced"):
        raise ValueError("projection must be 'unbalanced' or 'balanced'")
    if max_states <= 0:
        raise ValueError("max_states must be positive")

    o = n - v
    # Precompute b[oi] = C(v+oi, v) for fast vars computation
    b = [comb(v + oi, v) for oi in range(o + 1)]

    def vars_of(vec: Tuple[int, ...]) -> int:
        out = 1
        for oi in vec:
            out *= b[oi]
        return out

    # -----------------------------
    # Balanced projection
    # -----------------------------
    if projection == "balanced":
        evaluated = 0
        # We respect max_states as a budget even though balanced space has size (o+1).
        limit = min(o + 1, max_states)
        for oi in range(limit):
            evaluated += 1
            vec = [oi] * l
            res = compute_S_and_vars(n, v, m, l, vec, mode=mode)
            S = int(res["S"])
            vars_ = int(res["vars"])
            ok = (S > vars_) if strict else (S >= vars_)
            if ok:
                num_projected_blocks = (l if oi < o else 0)
                return {
                    "status": "found",
                    "mode": mode,
                    "projection": projection,
                    "strict": strict,
                    "o_kept_vec_sorted": vec,
                    "num_projected_blocks": num_projected_blocks,
                    "vars": vars_,
                    "S": S,
                    "rankpred": S,
                    "n_blocks": res["n_blocks"],
                    "rmax": res["rmax"],
                    "states_evaluated": evaluated,
                    "n_eqs": int(res["n_eqs"]),
                    "M": int(res["M"]),
                }

        return {
            "status": "not_found_within_budget",
            "mode": mode,
            "projection": projection,
            "strict": strict,
            "states_evaluated": evaluated,
            "max_states": max_states,
        }

    # --------------------------------
    # Unbalanced projection
    # --------------------------------
    start = tuple([0] * l)  # smallest vars configuration
    heap = [(vars_of(start), start)]
    seen = {start}

    evaluated = 0
    while heap and evaluated < max_states:
        vars_est, vec = heapq.heappop(heap)
        evaluated += 1

        res = compute_S_and_vars(n, v, m, l, list(vec), mode=mode)
        S = int(res["S"])
        vars_ = int(res["vars"])
        ok = (S > vars_) if strict else (S >= vars_)

        if ok:
            num_projected_blocks = sum(1 for oi in vec if oi < o)
            return {
                "status": "found",
                "mode": mode,
                "projection": projection,
                "strict": strict,
                "o_kept_vec_sorted": list(vec),
                "num_projected_blocks": num_projected_blocks,
                "vars": vars_,
                "S": S,
                "rankpred": S,
                "n_blocks": res["n_blocks"],
                "rmax": res["rmax"],
                "states_evaluated": evaluated,
                "n_eqs": int(res["n_eqs"]),
                "M": int(res["M"]),
            }

        # Expand: increase one component by 1 (then re-sort to keep canonical)
        for i in range(l):
            if vec[i] >= o:
                continue
            new = list(vec)
            new[i] += 1
            new.sort()
            newt = tuple(new)
            if newt in seen:
                continue
            seen.add(newt)
            heapq.heappush(heap, (vars_of(newt), newt))

    return {
        "status": "not_found_within_budget",
        "mode": mode,
        "projection": projection,
        "strict": strict,
        "states_evaluated": evaluated,
        "max_states": max_states,
    }


# ---------------------------------------------------------------------
# Full interface
# ---------------------------------------------------------------------

def best_projection_both(
    n: int,
    v: int,
    m: int,
    l: int,
    projection: str = "unbalanced",
    strict: bool = True,
    max_states: int = 20000,
) -> Dict[str, object]:
    """
    Run the optimization for both conservative and aggressive modes, under the
    chosen projection regime ("unbalanced" or "balanced").
    """
    cons = find_best_projection(
        n, v, m, l,
        mode="conservative",
        projection=projection,
        strict=strict,
        max_states=max_states,
    )
    aggr = find_best_projection(
        n, v, m, l,
        mode="aggressive",
        projection=projection,
        strict=strict,
        max_states=max_states,
    )
    return {"conservative": cons, "aggressive": aggr}


def _log2_binom(n: int, k: int) -> float:
    if not (0 <= k <= n):
        return float("-inf")
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2.0)


def cost_complexity(s) -> float:
    term1 = 2 * s  # Sparse matrix solver complexity
    term2 = 2 * math.log2(v+1)  # Density (v+1)^2
    q = 16     # log2(q^l)
    L = l * math.log2(q)
    field_cost = 3 * (2 * (L ** 2) + L)
    term3 = math.log2(field_cost)  # Field-arithmetic cost: 3 * (2*(log2(q^l))^2  +  log2(q^l))
    return round(term1 + term2 + term3, 0)

# ---------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------


if __name__ == "__main__":
    param_set_concrete = [
        [37, 17, 17, 2],
        [25,  8,  8, 3],
        [24,  5,  5, 4],
        [56, 25, 25, 2],
        [49, 11, 11, 3],
        [37,  8,  8, 4],
        [24,  5,  5, 5],
        [75, 33, 33, 2],
        [66, 15, 15, 3],
        [60, 10, 10, 4],
        [29,  6,  6, 5],
        [34,  4,  4, 4],
    ]

    # Search settings
    STRICT = True          # require S > vars (set False for S >= vars)
    MAX_STATES = 20000     # raise if "not_found_within_budget" occurs often

    # Projection regime:
    #   - "unbalanced" (default):multiset / sorted tuple search
    #   - "balanced": all-equal tuple search
    PROJECTION = "unbalanced"

    for v, o, m, l in param_set_concrete:
        n = v + o
        print(f"Case: v={v}, o={o}, m={m}, l={l}")

        cons = find_best_projection(
            n, v, m, l,
            mode="conservative",
            projection=PROJECTION,
            strict=STRICT,
            max_states=MAX_STATES,
        )
        aggr = find_best_projection(
            n, v, m, l,
            mode="aggressive",
            projection=PROJECTION,
            strict=STRICT,
            max_states=MAX_STATES,
        )

        print(f"[Conservative | {PROJECTION} projection]")
        print(f"  status               : {cons.get('status')}")
        if cons.get("status") == "found":
            print(f"  o_kept_vec_sorted     : {cons['o_kept_vec_sorted']}")
            print(f"  n_eqs        : {math.log(cons['n_eqs'], 2):0f}")
            print(f"  vars                  : {math.log(cons['vars'], 2):0f}")
            print(f"  rank         : {math.log(cons['S'], 2):0f}")
            print(f"  n_eqs - rank : {(math.log(cons['n_eqs'], 2) - math.log(cons['S'], 2)):0f}")
            print(f"  cost         : {cost_complexity(math.log(cons['vars'], 2))}")

        print(f"[Aggressive | {PROJECTION} projection]")
        print(f"  status               : {aggr.get('status')}")
        if aggr.get("status") == "found":
            print(f"  o_kept_vec_sorted     : {aggr['o_kept_vec_sorted']}")
            print(f"  n_eqs        : {math.log(aggr['n_eqs'], 2):0f}")
            print(f"  vars                  : {math.log(aggr['vars'], 2):0f}")
            print(f"  rank         : {math.log(aggr['S'], 2):0f}")
            print(f"  n_eqs - rank : {(math.log(aggr['n_eqs'], 2) - math.log(aggr['S'], 2)):0f}")
            print(f"  cost         : {cost_complexity(math.log(aggr['vars'], 2)):0f}")

        print()  # blank line between cases
