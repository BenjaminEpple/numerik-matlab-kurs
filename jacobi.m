% Input code

clear; clc;
format short;

% create A, a trilinear matrix
n = 5;
v = 10 * ones(n, 1);
uo = -1 * ones(n - 1, 1);
uu = 2 * ones(n - 1, 1);
D = diag(v, 0);
Do = diag(uo, 1);
Du = diag(uu, -1);
A = D + Do + Du;

% Rechte Seite 
b = ones(n, 1);

itermax = 10;

% calculate
x = do_jacobi(A, b, itermax);

% show result
disp(A);
disp(b);
disp(x);


function x = do_jacobi(A, b, itermax)
    
    % calculate Diagonal matrix D
    size = length(b);
    D = eye(size);
    for i=1:size
        D(i, i) = A(i, i);
    end
    
    % Zerlegung von A=D-L-R
    R = triu(A) - D;
    L = tril(A) - D;
    R = -1 * R;
    L = -1 * L;
    
    % Berechne Iterationsmatrix und Summanden d
    M = D^-1 * (L + R);
    d = D^-1 * b;
    
    % Startvektor
    x = ones(size, 1);
    
    for i=1:itermax
        x = iter_jacobi(M, d, x);
    end
end

    
function x = iter_jacobi(M, d, x_old)   
    x = M * x_old + d;
end 




