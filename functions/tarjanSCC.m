function scc = tarjanSCC(G)
    n = size(G, 1);
    index = zeros(n, 1); % Pour conserver l'index assigné à chaque nœud
    lowlink = zeros(n, 1); % Pour conserver le plus bas index atteignable
    onStack = false(n, 1); % Booléen pour vérifier si un nœud est sur la pile
    S = []; % Pile pour conserver les nœuds
    idx = 1; % L'index actuel
    scc = {}; % Pour conserver les composantes fortement connexes
    
    for i = 1:n
        if index(i) == 0
            [index, lowlink, onStack, S, idx, scc] = strongconnect(i, index, lowlink, onStack, S, idx, scc, G);
        end
    end
end

function [index, lowlink, onStack, S, idx, scc] = strongconnect(v, index, lowlink, onStack, S, idx, scc, G)
    index(v) = idx;
    lowlink(v) = idx;
    idx = idx + 1;
    S(end+1) = v; % Ajouter v à la pile S
    onStack(v) = true;
    
    for w = find(G(v,:)) % Pour chaque nœud w adjacent à v
        if index(w) == 0
            [index, lowlink, onStack, S, idx, scc] = strongconnect(w, index, lowlink, onStack, S, idx, scc, G);
            lowlink(v) = min(lowlink(v), lowlink(w));
        elseif onStack(w)
            lowlink(v) = min(lowlink(v), index(w));
        end
    end
    
    % Si v est une racine d'une composante fortement connexe
    if lowlink(v) == index(v)
        % Commencer une nouvelle composante fortement connexe
        i = length(S);
        component = [];
        while i > 0 && S(i) ~= v
            component(end+1) = S(i);
            onStack(S(i)) = false;
            i = i - 1;
        end
        component(end+1) = S(i);
        onStack(S(i)) = false;
        S(i:end) = []; % Retirer la composante de la pile
        scc{end+1} = component; % Ajouter la nouvelle composante
    end
end

