%% DDML demo

clc;
clear all;

addpath('activations');

%% LFW data: sparse sift (ssift)

load('data/lfw_sift.mat');

rand('state', 0);

dim_layer = [300 200 150];  % [300 200 150]   %[250 150 100] -0.8672   0.8672    0.0160 /
        
nL = length(dim_layer);     % number of layers

ux = ux(1:dim_layer(1),:);  %300 rows
n_sam = size(ux, 2);    % = 12000  ??????? 
x_mean = mean(ux, 2);     %????? ??????   500X1
ux = ux - repmat(x_mean , 1, n_sam);    %500*12000

h = fspecial('motion',31,1);
ux = imfilter(ux, h, 'replicate');

%% initializing
[Wo, bo] = initialize_Weights(dim_layer, 1);

%% ddml
matches = logical([pairs.match]); %1*6000
un = unique([pairs.fold]);  %10 only one value without repeat 

nfold = length(un); %10
t_acc = zeros(nfold, 1);% zero vector 10*1
for c = 1:nfold
    
    disp ('*************************************************************************************');
    disp ('*************************************************************************************');
    disp(['fold:  ' num2str(c)]);
    
    trainMask = [pairs.fold] ~= c; % 9 folds except fold ?= c
    testMask = [pairs.fold] == c;  % one fold only == c
    tr_idxa = idxa(trainMask); % ???? 
    tr_idxb = idxb(trainMask); % ????
    tr_matches = matches(trainMask); 
    
    [W, b ] = ddml_bp(tr_idxa, tr_idxb, tr_matches, ux, Wo, bo, nL);  
    ts_idxa = idxa(testMask);
    ts_idxb = idxb(testMask);
    ts_matches = matches(testMask);
    
    Xa = ux(:, ts_idxa);
    Xb = ux(:, ts_idxb);
    Xa = ddml_fp(Xa, W, b, nL);
    Xb = ddml_fp(Xb, W, b, nL);
    
    % cosine similarity
    sim = cos_sim(Xa, Xb);
    % accuracy
    [~, ~, acc] = ROCcurve(sim, ts_matches);
    t_acc(c) = acc;
end


disp([mean(t_acc) std(t_acc)]); % show mean accuracy
stderror= std( t_acc ) / sqrt( length( t_acc ));
disp ('SE= ');
disp (stderror);