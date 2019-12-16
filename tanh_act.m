function [a, grad] = tanh_act(X)

a = tanh(X);

if nargout > 1
    grad = (1-a) .* (1+a);
end

end