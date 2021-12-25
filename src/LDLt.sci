rand('seed', getdate("s"))

A = [6 15 55;15 55 225;55 225 979]
printf("\n# A:")
disp(A)

// 	2.4494897 	0.  		0.
// 	6.1237244 	4.1833001 	0.
// 	22.453656   20.916501   6.1101009

function [L] = llt(A)
	n = size(A)(1)
	D = eye(n, n)
	L = eye(n,n)
	if (n == size(A)(2))
		sum = 0
		for r = 1:n 
			for c = 1:r
				sum = 0
				for s = 1:(c-1)
					sum = sum + L(r,s) * L(c,s)
				end
				sum = A(r,c) - sum
				if (r == c)
					D(r,c) = sqrt(sum)
				else
					L(r,c) = sum / D(c,c)
				end
			end
		end
		L = L + D - eye(n ,n)
	else
		disp("The matrix A needs to be a square matrix")
	end
endfunction

function [L,D] = ldlt(A)
	n = size(A)(1)
	D = eye(n, n)
	L = eye(n,n)
	if (n == size(A)(2))
		for i = 1:n
			for j = 1:n
				sum = 0.0
				for k = 1:(j-1)
					sum = sum + L(j,k) * L(j,k) * D(k,k)
				end
				D(j,j) = A(j,j) - sum
				if i > j then
					sum = 0.0	
					for k = 1:(j-1)
						sum = sum + L(i,k) * L(j,k) * D(k,k)
					end
					L(i,j) = (A(i,j) - sum) / D(j,j)
				end
			end
		end
	else
		disp("The matrix A needs to be a square matrix")
	end
endfunction


[L] = llt(A)
printf("\n# L:")
disp(L)

printf("\n# Cholesky LLt:")
disp(L*L')

[L,D] = ldlt(A)
printf("\n# L:")
disp(L)
printf("\n# D:")
disp(D)

printf("\n# LDLt:")
disp(L*D*L')


printf("\n# A:")
disp(A)

printf("\n# difference A & LDLt:")
disp(norm(A) - norm(L*D*L'))

mean_cost_lu = 0.0
mean_cost_ldlt = 0.0
for s = 10:+2:50
	printf("\n-- Size %d: \t", s)
	timer()
	A = rand(s,s)
	[L,D] = ldlt(A)
	t = timer()
	printf("\tLDLt: %f sec", t)
	mean_cost_ldlt  = mean_cost_ldlt + t/s
	[L,U] = lu(A)
	printf("\tLU: %f sec\t--", t)
	mean_cost_lu = mean_cost_lu + t/s
end
printf("\nMean cost per size unit for:\n\tLDLt : %f seconds\n\tLU : %f seconds\n\n", mean_cost_ldlt, mean_cost_lu)
