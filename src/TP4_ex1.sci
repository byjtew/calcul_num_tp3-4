rand('seed', getdate("s"))

A = [6 15 55;15 55 225;55 225 979]
printf("\n# A:")
disp(A)

// 	2.4494897 	0.  		0.
// 	6.1237244 	4.1833001 	0.
// 	22.453656   20.916501   6.1101009

function [L, D] = llt(A)
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
	else
		disp("The matrix A needs to be a square matrix")
	end
endfunction


[L,D] = llt(A)
printf("\n# L:")
disp(L)
printf("\n# D:")
disp(D)

LD = L+D-eye(size(A)(1), size(A)(1))
printf("\n# Cholesky LLt:")
disp(LD*LD')
printf("\n# LDLt:")
disp(L*D*(L'))
printf("\n# A:")
disp(A)
