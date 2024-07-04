function sccGraphviz(G, grid_step, x_inf, y_inf, scc, filename)
n = size(G, 1);
w = sqrt(n);
dotStr = 'digraph G {\n';
% colors = [
%     "#ff0000",  % Red
%     "#00ff00",  % Green
%     "#0000ff",  % Blue
%     "#ffff00",  % Yellow
%     "#ff00ff",  % Magenta
%     "#00ffff",  % Cyan
%     "#808080",  % Gray
%     "#ff8000",  % Orange
%     "#800080",  % Purple
%     "#008080"   % Teal
%     ];
colors = [
    "#757474",  % light gray 808080
    ];
nodesAdded = {};

for i = 1:n
    for j = 1:n
        if G(i, j) == 1
            y = ceil(i/w);
            x = i - (y - 1) * w;
            y_prime = ceil(j/w);
            x_prime = j - (y_prime - 1) * w;
            nodeId_i = strrep(['_', num2str(x), '_', num2str(y)], '-', 'm');
            nodeId_j = strrep(['_', num2str(x_prime), '_', num2str(y_prime)], '-', 'm');
            [isIInSCC, ri] = findSCCIndex(i, scc);
            [isJInSCC, rj] = findSCCIndex(j, scc);
%             if isIInSCC && isJInSCC
                dotStr = [dotStr, nodeId_i, ' -> ', nodeId_j, ' [color="#0000001A"];\n'];
%             end
            
            if ~ismember(nodeId_i, nodesAdded)
                labelStr_i = [num2str(x), ',', num2str(y)];
                positionStr_i = [num2str(2*x), ',', num2str(2*y)];
                positionStr_i_float = [num2str(x_inf(1) + (x-1) * grid_step, '%.2f'), ',', num2str(y_inf(1) + (y-1) * grid_step, '%.2f')];
                if isIInSCC
                    dotStr = [dotStr, nodeId_i, ' [label="', positionStr_i_float, '", pos="', positionStr_i, '!", style=filled, fillcolor="', char(colors(mod(ri, length(colors)) +1)), '"];\n'];
                    nodesAdded{end+1} = nodeId_i;
                else
                 dotStr = [dotStr, nodeId_i, ' [label="', positionStr_i_float, '", pos="', positionStr_i, '!"];\n'];

                end
            end
            
            if ~ismember(nodeId_j, nodesAdded)
                labelStr_j = [num2str(x_prime), ',', num2str(y_prime)];
                positionStr_j = [num2str(2*x_prime), ',', num2str(2*y_prime)];
                positionStr_j_float = [num2str(x_inf(1) + (x_prime-1) * grid_step, '%.2f'), ',', num2str(y_inf(1) + (y_prime-1) * grid_step, '%.2f')];
                if isJInSCC
                    
                    dotStr = [dotStr, nodeId_j, ' [label="', positionStr_j_float, '", pos="', positionStr_j, '!", style=filled, fillcolor="', char(colors(mod(rj, length(colors)) +1)), '"];\n'];
                   
                    nodesAdded{end+1} = nodeId_j;
                else
                    dotStr = [dotStr, nodeId_j, ' [label="', positionStr_j_float, '", pos="', positionStr_j, '!"];\n'];
                end
                
            end
        end
    end
end

dotStr = [dotStr, '}'];
fileID = fopen(strcat('\\wsl$\Ubuntu\home\gwen\graphs\', filename, '.dot'), 'w');

fprintf(fileID, dotStr);
fclose(fileID);
[status, cmdout] = system(strcat('wsl neato -Tpdf /home/gwen/graphs/', filename, '.dot -o /home/gwen/graphs/', filename,'.pdf'));
[status, cmdout] = system(strcat('wsl neato -Tpng /home/gwen/graphs/', filename, '.dot -o /home/gwen/graphs/', filename,'.png'));
[status, cmdout] = system(strcat('wsl neato -Teps /home/gwen/graphs/', filename, '.dot -o /home/gwen/graphs/', filename,'.eps'));

end