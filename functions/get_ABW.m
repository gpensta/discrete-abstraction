function [A, B, W] = get_ABW(f, S, U)
    s_temp = S.interval(); x1 = s_temp(1); x2 = s_temp(2);
    u_temp = U.interval(); u1 = u_temp(1); u2 = u_temp(2);
    n = 2;
    m = 2;
    x0 = S.center;
    u0 = U.center;
    x_sym = sym('x', [n, 1], 'real');
    u_sym = sym('u', [m, 1], 'real');
    f_sym = f([x_sym; u_sym]);
    A_sym = jacobian(f_sym, x_sym);
    B_sym = jacobian(f_sym, u_sym);
    % Evaluate the Jacobian matrices at the equilibrium point
    A = double(subs(A_sym, [x_sym; u_sym], [x0; u0]));
    B = double(subs(B_sym, [x_sym; u_sym], [x0; u0]));
    diff = @(x) f([x(1); x(2); x(3); x(4)]) - (f([x0(1); x0(2); u0(1); u0(2)]) + A * ([x(1); x(2)] - x0) + B * ([x(3); x(4)] - u0));
    D = interval([x1.inf, x2.inf, u1.inf , u2.inf], [x1.sup, x2.sup, u1.sup , u2.sup]);
    tay = taylm(D,10,'x','linQuad');
    W = interval(diff(tay));
end