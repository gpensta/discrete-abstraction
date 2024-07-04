function x_final = apply_fn_iteratively(f, x, u_list)
    for i = 1:length(u_list)
        u = u_list(i);
        x = f(x, u);
    end
    x_final = x;
end