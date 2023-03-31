% example data

A = [1 5 4;
     -2 1 3;
     2 0 0];

b = [1 0 1]';
 
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
 
disp("End product of givens rotations:")
disp(A)

% compute solution x for Ax = b
x = A \ b;

disp("Result:")
disp(x)

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
 
 