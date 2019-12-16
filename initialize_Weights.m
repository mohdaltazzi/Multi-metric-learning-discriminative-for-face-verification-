function [Wo, bo] = initialize_Weights(dim_layer, mark)
%% initialize Weights: W, b.

nL = length(dim_layer) - 1;

Wo = cell(nL, 1);
bo = cell(nL, 1);

switch mark
    case 0
        % 0: zeros
        for k = 1:nL;
            Wo{k} = zeros(dim_layer(k+1), dim_layer(k));
            bo{k} = zeros(dim_layer(k+1), 1);
        end
    case 1
        % 1: eye
        for k = 1:nL;
            Wo{k} = eye(dim_layer(k+1), dim_layer(k));
            bo{k} = zeros(dim_layer(k+1), 1);
        end
    case 2
        % 2
        for k = 1:nL;
            r = sqrt(6) / sqrt(dim_layer(k+1) + dim_layer(k));
            Wo{k} = rand(dim_layer(k+1), dim_layer(k)) * 2 * r - r;
            bo{k} = zeros(dim_layer(k+1), 1);
        end       
    case 3
        % 3: [-0.5  0.5]
        for k = 1:nL;
            Wo{k} = 0.5 - rand(dim_layer(k+1), dim_layer(k));
            bo{k} = zeros(dim_layer(k+1), 1);
        end   
end

end