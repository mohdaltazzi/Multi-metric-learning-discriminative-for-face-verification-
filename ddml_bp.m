function [W, b ] = ddml_bp(idxa,idxb,matches, X, Wo, bo, nL)

%% DDML: back-propagation 

Th = 3;  % threshold: Th = 3 *




beta = 10;  % sharpness =10 *          


lambda = 1e-3;  % regularization: 1e-3 *





stepsize = 1e-4; % learning rate: 1e-4 *


           
maxiter = 50;      % max iteration: 50 *






lambda1 =    2.0612e-10;



obj = zeros(maxiter, 1);
mean_loss=zeros(maxiter, 1);
Xa = X(:, idxa);
Xb = X(:, idxb);
clear X;

[~, N] = size(Xa);

label = ones(N, 1);
label(~matches) = -1;

%% initializing

dW = cell(size(Wo));
db = cell(size(bo));
for i = 1:nL-1
    dW{i} = zeros(size(Wo{i}));
    db{i} = zeros(size(bo{i}));
end

W = Wo;
b = bo;
clear Wo bo;

%% update
abs_obj= zeros(1,50);
Obj_no = zeros(1,26);
v = zeros(1,N);
lossf = zeros(1,N);
lossff =zeros(1,N);
h1 = cell(nL, 1);
h2 = cell(nL, 1);
dh1 = cell(nL, 1);
dh2 = cell(nL, 1);
delta1 = cell(nL, 1);
delta2 = cell(nL, 1);
for iter = 1:maxiter
    idx = randperm(N);
    for i = 1:N   % N= 5400
        h1{1} = Xa(:, idx(i));
        h2{1} = Xb(:, idx(i));
        for j = 1:nL-1
            [h1{j+1}, dh1{j+1}] = tanh_act(W{j} * h1{j} + b{j});
            [h2{j+1}, dh2{j+1}] = tanh_act(W{j} * h2{j} + b{j});
           % [h1{j+1}, dh1{j+1}] = sigmoid_act(W{j} * h1{j} + b{j});
           % [h2{j+1}, dh2{j+1}] = sigmoid_act(W{j} * h2{j} + b{j});
            %[h1{j+1}, dh1{j+1}] = nonsaturate_sigmoid_act(W{j} * h1{j} + b{j});
            %[h2{j+1}, dh2{j+1}] = nonsaturate_sigmoid_act(W{j} * h2{j} + b{j});
            %[h1{j+1}, dh1{j+1}] = linear_act(W{j} * h1{j} + b{j});
            %[h2{j+1}, dh2{j+1}] = linear_act(W{j} * h2{j} + b{j});
        end
        
        dist = (h1{nL} - h2{nL})' * (h1{nL} - h2{nL});       
        x = 1 - label(idx(i)) * (Th - dist);  
        log_x = 1 / (1 + exp(-beta * x)); 
         
        
        lossf(idx(i))=log_x;
         
           if (log_x > lambda1 ) %------------------------------------>  V =0/1 if loss function < or > lambda
             v(idx(i)) = 0;
               log_x= log_x * v(idx(i));     
           else
             v(idx(i)) = 1;
              log_x= log_x * v(idx(i)) ;
           end 
           
        lossff(idx(i))=log_x;
        
        cur_lab = log_x * label(idx(i));
        delta1{nL} = cur_lab * ((h1{nL}-h2{nL}) .* (dh1{nL})); %Q10
        delta2{nL} = cur_lab * ((h2{nL}-h1{nL}) .* (dh2{nL})); %Q11
        for j = (nL-1):-1:2
            delta1{j} = (W{j}' * delta1{j+1}) .* (dh1{j});% Q12
            delta2{j} = (W{j}' * delta2{j+1}) .* (dh2{j}); % Q13
        end
        for j = 1:(nL-1)
           % dW{j} =  delta1{j+1} * h1{j}' + delta2{j+1} * h2{j}';
          % db{j} =   delta1{j+1} + delta2{j+1};
           
             dW{j} =  v(idx(i))*(delta1{j+1} * h1{j}' + delta2{j+1} * h2{j}');%Q8
             db{j} =  v(idx(i))*(delta1{j+1} + delta2{j+1}); %Q9
             W{j} = W{j} - stepsize * (dW{j} + lambda * W{j}); %Q16
             b{j} = b{j} - stepsize * (db{j} + lambda * b{j}); %Q17  
        end  
    end   % end N loops
    obj(iter) = com_obj(Xa, Xb, W, b, nL, label', lambda, Th, beta ,v ,lambda1);
     Obj_no(iter)= obj (iter);
     
    if (iter > 1) && (abs(obj(iter) - obj(iter-1)) < 1e-2) %<1e-2
      
        disp(iter);
        return;
    end
 
   lambda1 = lambda1 + 0.001 ;   % Increasing       
        
end
end

function obj = com_obj(Xa, Xb, W, b, nL, label, lambda, Th, beta,v,lambda1)
%% computing value of objective function

[~, N] = size(Xa);

t_h1 = ddml_fp(Xa, W, b, nL);
t_h2 = ddml_fp(Xb, W, b, nL);

diff = t_h1 - t_h2;
diff = diff .^ 2;
dist = sum(diff, 1);
x = 1 - label .* (Th - dist);
x = beta .* x;

J1 = sum(x(x > 10));
J2 = sum(log(1 + exp(x(x <= 10))));
J = (J1 + J2) / beta;
J = J / N;

r_W = 0;
r_b = 0;
for j = 1:(nL-1)
    r_t = norm(W{j}, 'fro');
    r_W = r_W + r_t * r_t ;
    r_t = norm(b{j}, 'fro');
    r_b = r_b + r_t * r_t;
end   
Spl = SPL_Req_Hard(v,lambda1);            %------------------> spl Hard req 

 obj = 0.5 * J + (0.5 * lambda) * (r_W + r_b) + Spl;
 


end