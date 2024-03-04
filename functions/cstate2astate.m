function astate = cstate2astate(cstate, xgrid, ygrid)
    % Extraire les coordonnées de cstate
    x = cstate(1);
    y = cstate(2);

    % Trouver les indices dans les grilles
    [i, j] = find(xgrid == x & ygrid == y);

    % Assurer qu'une correspondance unique est trouvée
    if length(i) ~= 1 || length(j) ~= 1
        error('Les coordonnées de cstate doivent correspondre à un point unique dans les grilles');
    end

    % Calculer astate
    w = size(xgrid, 2);
    astate = (i - 1) * w + j;
end