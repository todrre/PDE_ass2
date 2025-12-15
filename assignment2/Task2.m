m = 50;
[A, b] = get_system(m);

%compare solutions between gs and pcg and chol and \ using tic and tok
tol = 1e-6;
maxit = 100000;
x0 = zeros(size(b));

% Gather all times
times = struct();

tic
x_gs = gs(A, b, x0, tol, maxit);
times.gs = toc;

tic
x_pcg = pcg(A, b, tol, maxit);
times.pcg = toc;

tic
L = chol(A, 'lower');
x_chol = L' \ (L \ b);
times.chol = toc;

tic
x_backslash = A\b;
times.backslash = toc;

%disp('Solution using Gauss-Seidel method:');
%disp(x_gs);
disp(['Time taken by Gauss-Seidel: ', num2str(times.gs), ' seconds']);

%disp('Solution using Preconditioned Conjugate Gradient method:');
%disp(x_pcg);
disp(['Time taken by PCG: ', num2str(times.pcg), ' seconds']);

%disp('Solution using MATLAB Chol operator:');
%disp(x_chol);
disp(['Time taken by cholesky operator: ', num2str(times.chol), ' seconds']);


%disp('Solution using MATLAB backslash operator:');
%disp(x_backslash);
disp(['Time taken by backslash operator: ', num2str(times.backslash), ' seconds']);


%compare norms between solutions
norm_gs_pcg = norm(x_gs - x_pcg);
norm_gs_chol = norm(x_gs - x_chol);
norm_pcg_chol = norm(x_pcg - x_chol);
disp(['Norm between GS and PCG solutions: ', num2str(norm_gs_pcg)]);
disp(['Norm between GS and Chol solutions: ', num2str(norm_gs_chol)]);
disp(['Norm between PCG and Chol solutions: ', num2str(norm_pcg_chol)]);
