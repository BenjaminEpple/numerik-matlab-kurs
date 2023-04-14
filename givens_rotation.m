% example data

A = [1 5 4;
    -2 1 3;
    2 0 0];

b = [1 0 1]';


[A, b] = apply_givens(A, b);

disp("End product of givens rotations:")
disp(A)

% compute solution x for Ax = b
x = A \ b;

% compute solution x for Ax = b
x_own = solve_tri_matrix(A, b);

disp("Result:")
disp(x)

disp("Own Result:")
disp(x_own)





function [A, b] = apply_givens(A, b)

% Givens rotation
n = size(A);
for i = 1:n
    for j=i+1:n
        
        % determine s, c
        [s, c] = get_sc_vals(A, i, j);
        
        G = givens_matrix(n, i, j, c, s);
        A = G*A;
        b = G*b;
    end
end
end


function [s, c] = get_sc_vals(A, i, j)
r_ii = sqrt(A(i, i)^2 + A(j, i)^2);
s = A(j,i) / r_ii;
c = A(i,i) / r_ii;

end


function G = givens_matrix(n, i, j, c, s)

G = eye(n);

% set sine and cosine values
G(i, i) = c;
G(j, j) = c;
G(i, j) = s;
G(j, i) = -s;

end

function x = solve_tri_matrix(A, b)
n = size(A);
vector_size = n(1);
x = zeros(vector_size, 1);
for i=n:-1:1
    for j=n:-1:i
        denominator = A(i, i);
        sum = get_sum(A, b, i);
        x(i) = (b(i) - sum) / denominator;
    end
end

end


function sum = get_sum(A, x, i)
n = size(A);
sum = 0;
for j=n:-1:i+1
    sum = sum + A(i, j) * x(j);
end
end
