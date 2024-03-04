function [A, B] = linearize(f, x0, u0)
    % Define x and u as symbols
    n = length(x0);
    m = length(u0);
    
    x_sym = sym('x', [n, 1], 'real');
    u_sym = sym('u', [m, 1], 'real');
    
    % Convert the function handle to symbolic expression
    f_sym = f(x_sym, u_sym);

    % Calculate the Jacobian matrices A and B
    A_sym = jacobian(f_sym, x_sym);
    B_sym = jacobian(f_sym, u_sym);

    % Evaluate the Jacobian matrices at the equilibrium point
    A = double(subs(A_sym, [x_sym; u_sym], [x0; u0]));
    B = double(subs(B_sym, [x_sym; u_sym], [x0; u0]));

end
