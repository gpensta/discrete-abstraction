function x_final = apply_fn_iteratively(f, x, u_list)
    g = -9.81;
    l = 1;
    beta = 1; 
    m = 1; 
    dt = .1;
    for i = 1:length(u_list)
        u = u_list(i);
        x = f(x, u);
    end
    x_final = x;
end