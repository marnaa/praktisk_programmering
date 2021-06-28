Answer to the ppnm examproject 2 (Lanczos algorithm)
Here, in lanczos.c, we implement the Lanczos algorithm which is a way to represent a real and symmetric matrix. 
This representation conserves the high(est) eignevalues depending on the size of the representation,
the kyrlov supspace. 

The project relies on the eigenvalue methods (Jacobi diagonalisation) developed in the previous homework. 

The following tests were made to make sure the implementation was succesful. 

In out.matrix_test.txt we use a Krylov subspace of size rank(A)=10 and test that A=V T V^T
That V^T V = 1 and that T =V^T A V. We also test that T is tridiagonal, but since this is by construction, 
it is not surprising. 

Next up, still in out.matrix_test.txt we choose the kyrlov supspace to be of size 6 while Rank(A) = 10, still. Here we check T- V^T A V = 0 
and V^T V = 1 and again that T is tridiagonal. 

These show that the algorithm works as intended !!!. 

In out.eigval_check.txt we check that the eigenvalues of T indeed approximate the highest eigenvalues of A.
Therefore, increasing the size of the kyrlov subspace from 1 to 10 iteratively, we see how the approximation 
becomes better and better, until it represents A completely at a size of 10. Already pretty early we see that
the eigenvalue of T approximates the larges eigenvalue of A pretty well. 

The eigenvalues are found through the Jacobi diagonalisation developed in the eigenvalue homework. 

One can check through the added GS_arnoldi function, which includes all the previous vector in orthogonalisation,
that T (or H in the script) actually becomes tridiagonalized, and therefore that the GS_arnoldi function uses 
unecessary computation for the symmetric and real problem. 

All in all lanczos algorithm is succesfully implemented and aproximates the highest eigenvalues af A on a smaller
subspace. (Self-evaluation [10/10])
  