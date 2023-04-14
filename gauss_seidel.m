% Input code

format short;

% create A, a trilinear matrix
n = 5;

% Vektoren für das Befüllen der Tridiagonalmatrix
v = 10 * ones(n, 1);
uo = -1 * ones(n - 1, 1);
uu = 2 * ones(n - 1, 1);

D = diag(v, 0);
D_oben = diag(uo, 1);
D_unten = diag(uu, -1);

A = D + D_oben + D_unten;

% Rechte Seite 
b = ones(n, 1);

itermax = 10;

% calculate
x = do_gauss_seidel(A, b, itermax);

% show result
disp(A);
disp(b);
disp(x);

function x = do_gauss_seidel(A, b, itermax)
    
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
    
    d = D^-1 * b;
    
    % Startvektor
    x = ones(size, 1);
    
    for i=1:itermax
        x = iter_gauss_seidel(D, L, R, d, x);
    end
end


function x = iter_gauss_seidel(D, L, R, d, x_old)
    
    % calculate needed matrices
    D_inverse_times_L =  D^-1 * L;
    D_inverse_times_R =  D^-1 * R;
    
    % result vector
    x = zeros(length(x_old), 1);
    
    for i=1:length(x)
        a = D_inverse_times_L(i,:) * x;
        b = D_inverse_times_R(i,:) * x_old;
        x(i) = d(i) + a + b;
    end
end 



