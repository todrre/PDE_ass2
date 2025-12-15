function x = gs(A, b, x0, tol, maxit)
%GS Solve the linear system Ax = b using the Gauss-Seidel method
% Inputs:
%   A     - Coefficient matrix (n x n)
%   b     - Right-hand side vector (n x 1)
%   x0    - Initial guess for the solution (n x 1)
%   tol   - Tolerance for convergence (scalar)
%   maxit - Maximum number of iterations (integer)
% Output:
%   x     - Approximate solution vector (n x 1)

% compare norm_2(Ax- b)/norm_2(b) < tol
% x^(k+1) = (L + D)^(-1) (b - Ux^(k))

% Decompose A into L, D, U
L = tril(A, -1);
D = diag(diag(A));
U = triu(A, 1);

x = x0;
for k = 1:maxit
    x_old = x;
    x = (L + D) \ (b - U * x_old);

    % Check for convergence
    if norm(A*x - b, 2) / norm(b, 2) < tol
        return;
    end
end

fprintf('Gauss-Seidel did not converge within the maximum number of iterations (%d).\n', maxit);
end

