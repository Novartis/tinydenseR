# Numerical Stability Test Plan for tinydenseR

## Executive Summary

**Objective**: Add lightweight, principled unit tests to ensure tinydenseR produces reliable and reproducible results across varying computational environments and input perturbations.

**Scope**: 4 test categories covering ~150 lines of test code in `test-numerical-stability.R`

**Time estimate**: 4 hours implementation + review

**Regulatory context**: Critical for Novartis validation requirements and publication reproducibility

---

## Test Categories

### 1. Seed Reproducibility ✓ CRITICAL
**Goal**: Verify bit-identical results with fixed `set.seed()`

**Tests**:
- Full pipeline determinism (setup → landmarks → graph → map)
- Clustering determinism with same seed
- Verify UMAP embedding, cluster IDs, and fdens matrices are identical

**Why**: Regulatory compliance requires reproducible statistical results. Non-determinism breaks trust.

**Coverage**: ~40 lines

---

### 2. Input Invariance ✓ CRITICAL
**Goal**: Ensure robustness to benign transformations

**Tests**:
- **Perturbation robustness**: Add noise (σ = 1e-8) → verify >90% Jaccard overlap in landmarks
- **Permutation invariance**: Shuffle cell order → verify identical fdens aggregates
- **Weight validity**: Assert fuzzy graph weights are non-negative, finite, sum correctly

**Why**: Small floating-point differences (BLAS implementations, rounding) shouldn't flip DA calls.

**Coverage**: ~60 lines

---

### 3. Edge Cases ✓ IMPORTANT
**Goal**: Graceful handling of boundary conditions

**Tests**:
- Minimum viable dimensions (n=2 samples, k=2 markers)
- Straggler absorption (clusters < threshold are merged)
- Empty/isolated nodes in kNN graph
- Extreme marker values (near-zero, very large counts)

**Why**: Real datasets have outliers, small batches, and technical artifacts. Code should fail gracefully or handle correctly.

**Coverage**: ~40 lines

---

### 4. CI Robustness ✓ NICE-TO-HAVE
**Goal**: Verify stable behavior across R versions and platforms

**Tests**:
- RNG consistency (different RNGkind() settings)
- Sparse matrix numerical properties (symmetry, value ranges)
- Cross-platform sanity checks (tagged with `skip_on_cran()` if slow)

**Why**: GitHub Actions runs multiple R versions; CRAN tests on Solaris/Windows. Avoid platform-specific failures.

**Coverage**: ~20 lines

---

## Implementation Strategy

### Files Modified
1. **`tests/testthat/test-numerical-stability.R`** (NEW)
   - All 4 test categories
   - Well-commented with rationale
   - Uses existing `helper-test-utils.R` functions

2. **`tests/testthat/helper-test-utils.R`** (minor additions if needed)
   - Add fixture for perturbed data generation
   - Add Jaccard similarity helper

### Test Execution
- All tests run in **<10 seconds** on typical hardware
- Expensive tests use `skip_on_cran()` to avoid CRAN timeout
- GitHub Actions CI runs full suite on each PR

### Coverage Goals
- Aim for **>80% line coverage** in core numerical modules:
  - `R/lm.graph.embed.R` (graph construction, clustering, mapping)
  - `R/landmarks.R` (landmark selection, PCA)
  - `R/stats.model.R` (density aggregation, DA testing)

---

## What This Plan DOES NOT Include

❌ **Rmpfr/arbitrary precision tests** - Overkill; double precision is sufficient  
❌ **Sub-machine-epsilon checks** - Real data never has distances <1e-12  
❌ **Design matrix rank testing** - That's `lm()`/`limma`'s job, not ours  
❌ **Testing external dependencies** - Trust `uwot`, `igraph`, `limma` to handle their numerics  
❌ **Exhaustive softmax stabilization** - `uwot` already uses stabilized UMAP; verify outputs, don't reimplement  

---

## Acceptance Criteria

**Definition of Done**:
1. ✅ All 4 test categories implemented and passing
2. ✅ Tests run in <10 seconds locally
3. ✅ GitHub Actions CI passes on all R versions
4. ✅ Code coverage report shows >80% for numerical modules
5. ✅ Documentation includes brief comment on reproducibility guarantees

**Validation**:
- Run `devtools::test()` locally → all pass
- Run `devtools::check()` → no new NOTEs/WARNINGs
- Trigger CI on PR → green checkmarks
- Review Codecov diff → acceptable coverage increase

---

## Open Questions for Review

1. **Tolerance levels**: Is `1e-12` appropriate for fdens equality, or should we use `1e-10`?
2. **Perturbation noise**: Is σ = 1e-8 realistic, or should we test σ = 1e-6 (larger but still tiny)?
3. **Jaccard threshold**: Is 90% overlap too strict for landmark stability? Should we use 85%?
4. **CI platforms**: Do we need to test on Windows/Mac explicitly, or is Linux sufficient?

---

## References

- **Test file**: `tests/testthat/test-numerical-stability.R` (attached)
- **Helper functions**: `tests/testthat/helper-test-utils.R` (existing)
- **Numerical computing best practices**: Higham, "Accuracy and Stability of Numerical Algorithms"
- **R package testing**: Wickham & Bryan, "R Packages" (2nd ed.), Chapter 13

---

## Next Steps

1. **Review this plan** with team (15 min)
2. **Implement tests** (~4 hours)
3. **Run locally** and fix any failures
4. **Create PR** with coverage report
5. **Merge** after approval

**Assignee**: [TBD]  
**Target completion**: [TBD]  
**Priority**: High (blocks v1.0 release)
