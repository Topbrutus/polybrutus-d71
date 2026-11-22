# Polybrutus 71D: Exact 76-Point Euclidean Metric and Sixty-Fivefold Ideal Graph Mode
**Human scientific author:** Gabriel St-Pierre, independent researcher.
**Material AI assistance:** OpenAI ChatGPT (Astra) assisted the exact programs, verification, manuscript, metadata and publication; AI is not the scientific author.

## 1. Conservative Euclidean extension
Embed the 75 unit vectors of the immutable D70 metric as q_i=(p_i,0), and append one orthogonal unit vector q_76. The 2775 inherited labelled squared distances remain unchanged; the 75 new distances are all 2R². The exhaustive total is 2850 unordered labelled pairs. At R=1, use the original seven-point canonical Gram block (rank 3), a regular five-point centered simplex (rank 4, off-diagonal -1/4), and 64 independent orthogonal unit axes (rank 64). The orthogonal sum has rank 71 and is positive semidefinite. The canonical antipodal pair ensures that the origin belongs to the affine hull, so the minimal affine dimension is 71, with 72 storage coordinates in this representation.

## 2. Exact metric invariants
With G the rational unit-diagonal Gram, D=2(11ᵀ-G), H=I-11ᵀ/76, the centered Gram B=-HDH/2=HGH is positive semidefinite and has rank 71. The labelled sum of squared distances is 22847/4. The exact histogram k=4d²/R² is 4:6, 7:2, 8:2819, 9:3, 10:10, 12:6, 13:2, 16:2 (2850 total).

## 3. Separately stipulated ideal graph
For equal scalar unit masses and complete-graph weights w_ij=1/d_ij² at R=1, the added point connects to every previous vertex with weight 1/2. The old Laplacian block is shifted by I/2. On old zero-sum modes eigenvalues shift by 1/2; the added contrast mode is 38 and the constant mode remains zero. Thus the characteristic polynomial satisfies:
  chi76(z)=z*(z-38)*chi75(z-1/2)/(z-1/2)
  chi76(z)=z*(z-38)^65*(z-75/2)^4*f6(z-69/2)
where f6(z)=chiP7(z)/z. Rational rank-nullity verifies 38 with multiplicity 65 and 75/2 with multiplicity 4. These are dimensionless model eigenvalues, not measured hertz.

## 4. Reproducibility and falsification
Two separately implemented exact rational programs construct the Gram or labelled distances independently, compare every one of 2850 labelled pairs and all 2775 inherited pairs, and compute centered Gram rank and graph nullities. Four adversarial controls reject corrupted distances, a duplicate final point, a tilted last point and histogram-preserving relabelling. The inherited D70 ZIP CRC, SHA256 and original certificate bytes are checked in both programs.

## 5. Immutable provenance and scope
Direct parent: Topbrutus/polybrutus-d70@v0.1.0, commit 0c56da5cd75613453226fc8ed252dacc01bab22b, DOI 10.5281/zenodo.24433371. The complete original scientific ZIP and JSON certificate are embedded unchanged; older nested originals remain untouched. An inaccurate D46 Git commit label in the historical D47 certificate was disclosed since D48; no historical original is rewritten. This synthetic metric and separately chosen ideal graph do not establish physical spacetime dimensions, observed frequencies, new particles, a novel universal theorem or independent human peer review.

## References
St-Pierre, Gabriel (2026). Polybrutus 70D: Exact 75-Point Euclidean Metric and Sixty-Fourfold Ideal Graph Mode. Zenodo. DOI: 10.5281/zenodo.24433371. https://doi.org/10.5281/zenodo.24433371

St-Pierre, Gabriel (2026). Polybrutus 7D: An Exact Twelve-Point Orthogonal Completion, the 12-35-37 Identity, and an Oscillator Model. Zenodo. DOI: 10.5281/zenodo.24407952. https://doi.org/10.5281/zenodo.24407952
