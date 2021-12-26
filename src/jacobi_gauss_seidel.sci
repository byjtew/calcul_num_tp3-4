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

function [A] = poisson_1d_dense(n) 
    A = zeros(n,n)
    A(1,1) = 2
    for i = 2:n
        A(i-1, i) = -1
        A(i, i-1) = -1
        A(i,i) = 2
    end
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

function [x, nb_iter] = jacobi(A, b)
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
    printf("<> jacobi(A[%d,%d], b) in %d iterations.\n", n, n, nb_iter)
endfunction

function [x, nb_iter] = gauss_seidel(A, b)
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
    printf("<> gauss_seidel(A[%d,%d], b) in %d iterations.\n", n, n, nb_iter)
endfunction

function [rho] =  spectral_radius(A_dense)
    rho = max(real(spec(A_dense)))
endfunction

function [mini, maxi] =  min_max_eigenvalue(A_dense)
    eigenvalues = real(spec(A_dense))
    maxi = max(eigenvalues)
    mini = min(eigenvalues)
endfunction

function [x, nb_iter] = richardson(A, b, alpha)
    nb_iter = 0
    x_exact = A\b
    n = size(A_dense)(1)
    x = zeros(n, 1)
    D_p_L = tril(A)
    U = triu(A, +1)

    max_iter = 5000

    // alpha_optimal = .5 * (lambda_min(A) + lambda_max(A))
    
    while abs(norm(A\b) - norm(x)) >= 1e-7 && nb_iter < max_iter then
        x = x + alpha * (b - A * x)
        //printf("\nDifference is: %f ", abs(norm(x_exact) - norm(x)))
        nb_iter = nb_iter + 1
    end
    printf("<> richardson(A[%d,%d], b, alpha=%f) in %d iterations.\n", n, n, alpha, nb_iter)
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


[x, nb_iter_jacobi] = jacobi(A_dense, b)
printf("\n# x_jacobi:")
disp(x)

[x, nb_iter_gauss_seidel] = gauss_seidel(A_dense, b)
printf("\n# x_gauss_seidel:")
disp(x)

[x, nb_iter_richardson] = richardson(A_dense, b, .5)
printf("\n# x_richardson:")
disp(x)

x_exact = A_dense\b
printf("\n# x_exact:")
disp(x_exact)

eigenvalues = spec(A_dense)
printf("\n# eigenvalues:")
disp(eigenvalues)

fp = mopen("jacobi_gauss_seidel_richardson_results.dat", "wb")
nb_repetitions = 10
for n = 3:25
    mean_jacobi = 0.0
    mean_gauss_seidel = 0.0
    mean_richardson = 0.0
    for r = 1:nb_repetitions
        //A_dense = tri_diag_to_dense(random_strictly_dominant_tri_diag(n))
        A_dense = poisson_1d_dense(n)
        b = rand(n,1)

        [x, nb_iter_jacobi] = jacobi(A_dense, b)
        mean_jacobi = mean_jacobi + nb_iter_jacobi

        [x, nb_iter_gauss_seidel] = gauss_seidel(A_dense, b)
        mean_gauss_seidel = mean_gauss_seidel + nb_iter_gauss_seidel

        [min_eigen, max_eigen] = min_max_eigenvalue(A_dense)
        alpha = 2.0/(max_eigen+min_eigen)
        [x, nb_iter_richardson] = richardson(A_dense, b, alpha)        
        mean_richardson = mean_richardson + nb_iter_richardson
    end
    printf("%f %%\n", 100.0*n/20.0)
    mfprintf(fp, "%d %f %f %f\n", n, mean_jacobi/nb_repetitions, mean_gauss_seidel/nb_repetitions, mean_richardson/nb_repetitions)
end
mclose(fp)