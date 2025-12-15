A = [4 1 0 0;
     1 4 1 0;
     0 1 4 1;
     0 0 1 4];

b = [1; 1; 1; 1];
tol = 1e-6;
maxit = 100;
x0 = zeros(size(b));
x_gs = gs(A, b, x0, tol, maxit);

x = A\b;

disp('Solution using Gauss-Seidel method:');
disp(x_gs);

disp('Solution using MATLAB backslash operator:');
disp(x);