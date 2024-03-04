function [isInSCC, sccIndex] = findSCCIndex(i, scc)
    isInSCC = false;
    sccIndex = -1;
    
    % Parcourir chaque composante dans scc
    for idx = 1:length(scc)
        % Vérifier si i est dans la composante actuelle
        if ismember(i, scc{idx})
            isInSCC = true;
            sccIndex = idx;
            break; % Sortir de la boucle une fois i trouvé
        end
    end
end
