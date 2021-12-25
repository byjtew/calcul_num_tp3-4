rand('seed', getdate("s"))

function [A_dense] = tri_diag_to_dense(A_tri_diag)
	n = size(A_tri_diag)(2)
	A_dense = zeros(n,n)
	for i = 1:n-1
		A_dense(i+1,i) = A_tri_diag(3,i)
		A_dense(i,i) = A_tri_diag(2,i)
		A_dense(i,i+1) = A_tri_diag(1,i+1)
	end
	A_dense(n,n) = A_tri_diag(2,n)
endfunction

function [L,U] = split_dense_in_LU(A_dense)
	n = size(A_dense)(1)
	L = eye(n,n)
	U = zeros(n,n)

	for i = 1:n-1
		L(i+1,i) = A_dense(i+1,i)
		U(i,i) = A_dense(i,i)
		U(i,i+1) = A_dense(i,i+1)
	end
	U(n,n) = A_dense(n,n)
endfunction

function [A] = random_tri_diag_matrix(n)
	A = rand(3,n)
	A(1) = 0
	A(3*n) = 0
endfunction


n=6
A = random_tri_diag_matrix(n)
printf("\n# A:")
disp(A)

A_dense = tri_diag_to_dense(A)


function [retA] = lu_tri_diag(A_tri_diag)
	n = size(A_tri_diag)(2)
	retA = A_tri_diag
	if (3 == size(A)(1))
		for i = 2:n
			retA(3,i-1) = retA(3,i-1) / retA(2,i-1)
			retA(2,i) = retA(2,i) - retA(3,i-1) * retA(1,i)
		end
	else
		disp("The matrix A needs to be in tri-diagonal format")
	end
endfunction


[LU_A] = lu_tri_diag(A)
printf("\n# LU_A:")
disp(LU_A)

[L,U] = split_dense_in_LU(tri_diag_to_dense(LU_A))

printf("\n# L:")
disp(L)
printf("\n# U:")
disp(U)

printf("\n# L*U:")
disp(L*U)

printf("\n# A_dense:")
disp(A_dense)

printf("\n# difference A_dense & L*U:")
disp(norm(A_dense) - norm(L*U))

mean_cost = 0.0
for s = 250:+250:10000
	printf("\n-- Size %d: ", s)
	timer()
	[LU_A] = lu_tri_diag(random_tri_diag_matrix(s))
	t = timer()
	mean_cost = mean_cost + t/s
	printf(" %f seconds --", t)
end
printf("\nmean cost per size unit: %f seconds\n\n", mean_cost)
