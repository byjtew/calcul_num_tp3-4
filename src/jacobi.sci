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

function [A] = random_tri_diag(n)
	A = rand(3,n)
	A(1) = 0
	A(3*n) = 0
endfunction


function [A] = random_strictly_dominant_tri_diag(n)
	A = rand(3,n)
    A(2, 1) = 2.0 * A(1, 2)
    A(1, 1) = 0
	A(3, n) = 0
    for i = 2:(n-1)
         A(2, i) = 2.0 * (abs(A(1, i+1)) + abs(A(3, i-1)))
    end
    A(2, n) = 2.0 * A(3, n-1)
	
endfunction

function [A] = random_strictly_dominant_dense(n)
	A = rand(n,n)
	for i = 1:n
        sum = 0.0
        for j = 1:(i-1)
            sum = sum + abs(A(i, j))
        end
        for j = (i+1):n
            sum = sum + abs(A(i, j))
        end
        A(i, i) = 2.0*sum
    end

endfunction

function [x] = jacobi(A, b)
    printf("\n> jacobi(...)")
    nb_iter = 0
    x_exact = A\b
    n = size(A_dense)(1)
    x = zeros(n, 1)
    D = eye(n, n)
    for i  = 1:n
        D(i, i) = A(i, i)
    end
    E = -tril(A, -1)
    F = -triu(A, +1)
    
    while abs(norm(A\b) - norm(x)) >= 1e-7 then
        x = inv(D) * ((E+F) * x + b)
        //printf("\nDifference is: %f ", abs(norm(A\b) - norm(x)))
        nb_iter = nb_iter + 1
    end
    printf("\n< jacobi(...) in %d iterations.\n", nb_iter)
endfunction

function [x] = gauss_seidel(A, b)
    printf("\n> gauss_seidel(...)")
    nb_iter = 0
    x_exact = A\b
    n = size(A_dense)(1)
    x = zeros(n, 1)
    D_p_L = tril(A)
    U = triu(A, +1)
    
    while abs(norm(A\b) - norm(x)) >= 1e-7 then
        x = inv(D_p_L) * (b - U * x)
        //printf("\nDifference is: %f ", abs(norm(x_exact) - norm(x)))
        nb_iter = nb_iter + 1
    end
    printf("\n< gauss_seidel(...) in %d iterations.\n", nb_iter)
endfunction

function [x] = richardson(A, b, alpha)
    printf("\n> richardson(.., alpha=%f)", alpha)
    nb_iter = 0
    x_exact = A\b
    n = size(A_dense)(1)
    x = zeros(n, 1)
    D_p_L = tril(A)
    U = triu(A, +1)

    // alpha_optimal = .5 * (lambda_min(A) + lambda_max(A))
    
    while abs(norm(A\b) - norm(x)) >= 1e-7 then
        x = x + alpha * (b - A * x)
        //printf("\nDifference is: %f ", abs(norm(x_exact) - norm(x)))
        nb_iter = nb_iter + 1
    end
    printf("\n< richardson(.., alpha=%f) in %d iterations.\n", alpha, nb_iter)
endfunction


n=6
A = random_strictly_dominant_tri_diag(n)
disp(A)
A_dense = tri_diag_to_dense(A)

printf("\n# A:")
disp(A_dense)

b = rand(n,1)
printf("\n# b:")
disp(b)


[x] = jacobi(A_dense, b)
printf("\n# x_jacobi:")
disp(x)

[x] = gauss_seidel(A_dense, b)
printf("\n# x_gauss_seidel:")
disp(x)

[x] = richardson(A_dense, b, .5)
printf("\n# x_richardson:")
disp(x)

x_exact = A_dense\b
printf("\n# x_exact:")
disp(x_exact)