function X = ddml_fp(X, W, b, nL)
%% DDML forward ropagation 

N = size(X, 2);
for j = 1:nL-1
    X = tanh_act(W{j} * X + b{j}*ones(1, N));
    %X = sigmoid_act(W{j} * X + b{j}*ones(1, N));
   %X =  nonsaturate_sigmoid_act(W{j} * X + b{j}*ones(1, N));
  
   % X =  linear_act(W{j} * X + b{j}*ones(1, N));
end