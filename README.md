# QR Decomposition

This repository holds code pieces that implements QR factorization routines.

Such routines decompose some matrix $A \in \mathbb{R}^{m \times n}$ into an orthonormal matrix $Q$ and an upper-triangular matrix $R$ such that $A = QR$.

## Implemented Algorithms

### Gram-Schmidt orthogonalization
For each column of $A$, form a new basis by subtracting the column vector's projection onto all previously found bases.

### Givens rotation
Rotate two rows per inner iteration to eliminate non-zero elements below the diagonal in $R$. The implementation uses in-place modification of Givens rotation matrix $G$ for better performance.

### Householder transformation
Find the hyperplane that can induce a reflected vector of column vectors in $A$ that has no non-zero element below the diagonal in $R$.
