function x_final = apply_fn_iteratively(x, u_list)
    g = -9.81;
    l = 1;
    beta = 1; 
    m = 1; 
    dt = .1;
    for i = 1:length(u_list)
        u = u_list(i);
        x = [x(1) + dt * x(2); 
             x(2) + dt * (-g/l*sin(x(1)) - beta/(m*l^2)*x(2) + u)];
    end
    x_final = x;
end