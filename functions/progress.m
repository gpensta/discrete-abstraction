function progress(i)
    if i > 1
        fprintf(repmat('\b', 1, length(num2str(i-1)))); % Supprimer les caractères précédents
    end
    % Afficher le nouvel indice
    fprintf('%d', i);
end