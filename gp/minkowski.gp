/* minkowski.gp: reduction to the fundamental domain of \H_2/\Gamma_2*/

/* Minkowski-reduce a symmetric real positive definite 2x2 matrix:
   find U, M' s.t. M' = U M tU is Minkowski-reduced. */
congr(M, U) = U * M * U~;

minkowski_reduced(M) =
{
   M[2,2] >= M[1,1] && M[1,1] >= 2 * M[1,2] && M[1,2] >= 0;
}  

minkowski(M) =
{
   my(U, t, N);
   t = 1; \\true
   U = [1, 0; 0, 1];
   while (t,
	  if (2 * abs(M[1,2]) <= abs(M[1,1]),
	      if (abs(M[1,1]) <= abs(M[2,2]),
		  if (M[1,2] <= 0,
		      N = [1, 0; 0, -1];
		      M = congr(M, N);
		      U = N * U;
		     );
		  t = 0; \\end while loop
		  , \\else: |M[1,1]| > |M[2,2]|
		  N = [0, 1; -1, 0];
		  M = congr(M, N);
		  U = N * U;
		 );
	     );
	  N = [1, 0; -round(M[1,2]/M[1,1]), 1];
	  M = congr(M, N);
	  U = N * U;
	 );
   [U, M];
}

/* Handling sp4 matrices */
id(n) = matrix(n) + 1;
zero(n) = matrix(n);
sp4_from_blocks(a, b, c, d) =
{
   my(ac, bd);
   ac = matconcat([a~, c~])~;
   bd = matconcat([b~, d~])~;
   matconcat([ac, bd]);
}
blocks_from_sp4(gma) = [gma[1..2,1..2], gma[1..2,3..4], gma[3..4,1..2], gma[3..4,3..4]];
check_in_sp4(gma) =
{
   my(a, b, c, d);
   [a, b, c, d] = blocks_from_sp4(gma);
   a~ * c == c~ * a && b~ * d == d~ * b && a~ * d - c~ * b == matid(2);
}
sp4_action(gma, tau) =
{
   my(a, b, c, d);
   [a, b, c, d] = blocks_from_sp4(gma);
   /* if (!check_in_sp4(gma),
       error(gma, "is not in Sp_4(ZZ)");
      ); */
   (a * tau + b) * (c * tau + d)^(-1);
}

/* Compute \gamma and \tau' s.t. tau'\in\F_2 and \tau' = \gamma\tau */
reduce_to_F2(tau) =
{
   my(gma, U, N, j, a, b, c, d, FORGET);
   gma = matid(4);
   t = 1; \\true
   while (t,
	  \\Minkowski-reduce imaginary part
	  [U, FORGET] = minkowski(imag(tau));
	  N = sp4_from_blocks(U, zero(2), zero(2), (U^(-1))~);
	  tau = sp4_action(N, tau);
	  gma = N * gma;
	  \\reduce real part
	  [N, tau] = reduce_real_part(tau);
	  gma = N * gma;
	  \\end loop ?
	  t = 0;
	  \\check condition for the 19 test matrices
	  for (j=1, 19,
	       N = F2_test_matrix(j);
	       [a, b, c, d] = blocks_from_sp4(N);
	       if (abs(matdet(c * tau + d)) < 1,
		   t = 1; \\loop again
		   tau = sp4_action(N, tau);
		   gma = N * gma;
		  );
	      );
	 );
   [gma, tau];
}

/* Compute \gamma and \tau' such that \tau' = \gamma\tau has
   small real part */
reduce_real_part(tau) =
{
   my(gma);
   gma = sp4_from_blocks(matid(2), -round(real(tau)), zero(2), matid(2));
   [gma, sp4_action(gma, tau)];
}

/* F2 test matrices */
F2_test_matrix(j) = globalvar_F2_test_matrices[j];
globalvar_F2_test_matrices =
{
   [sp4_from_blocks(zero(2), -matid(2), matid(2), zero(2)),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [1, 0; 0, 0]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [-1, 0; 0, 0]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [0, 0; 0, 1]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [0, 0; 0, -1]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), matid(2)),
    sp4_from_blocks(zero(2), -matid(2), matid(2), -matid(2)),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [-1, 0; 0, 1]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [1, 0; 0, -1]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [0, 1; 1, 0]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [0, -1; -1, 0]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [1, 1; 1, 0]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [-1, -1; -1, 0]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [0, 1; 1, 1]),
    sp4_from_blocks(zero(2), -matid(2), matid(2), [0, -1; -1, -1]),
    sp4_from_blocks(matid(2), -matid(2), [1, 0; 0, 0], [0, 0; 0, 1]),
    sp4_from_blocks(matid(2), -matid(2), [0, 0; 0, 1], [1, 0; 0, 0]),
    sp4_from_blocks(matid(2), zero(2), [1, -1; -1, 1], id(2)),
    sp4_from_blocks(-matid(2), zero(2), [1, -1; -1, 1], -id(2))
   ]
}

F2_reduced(tau) =
{
   my(b_real, b_minkowski, b_test_matrices, j, a, b, c, d);
   b_real = abs(real(tau[1,1])) <= 1/2 && abs(real(tau[1,2])) <= 1/2
	    && abs(real(tau[2,2])) <= 1/2;
   b_minkowski = minkowski_reduced(imag(tau));
   b_test_matrices = 1;
   for (j=1, 19,
	[a, b, c, d] = blocks_from_sp4(F2_test_matrix(j));
	if (abs(matdet(c * tau + d)) < 1,
	    b_test_matrices = 0;
	   );
       );
   b_real && b_minkowski && b_test_matrices;
}

/* given tau1 and tau2 that we suppose equal in H_2/sp4, find
   gamma such that tau2 = gamma tau1.
   WARNING: this may fail if computations are inexact and tau lives
   on the boundary of the fundamental domain. */
find_sp4(tau1, tau2) =
{
   my(gma1, gma2, FORGET);
   [gma1, FORGET] = reduce_to_F2(tau1);
   [gma2, FORGET] = reduce_to_F2(tau2);
   gma2^(-1) * gma1;
}
