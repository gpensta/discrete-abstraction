function bool = existsTransition(A, B, alpha, x0, target, Uref)
    max_u1 = interval(Uref).sup(1);
    max_u2 = interval(Uref).sup(2);
    u = sdpvar(2, 1);
    Constraints = [u(1) <= max_u1, u(1) >= -max_u1, 
                   u(2) <= max_u2, u(2) >= -max_u2, 
                   norm(A * x0 + B * u - target, 'inf') <= alpha];
    options = sdpsettings('verbose', 0);
    diagnostics = optimize(Constraints, [], options);
    if diagnostics.problem == 0
        bool = 1;
    else 
        bool = 0;
    end
end
