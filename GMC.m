%
%% min sum_v{sum_i{||x_i - x_j||^2*s_ij + alpha*||s_i||^2} + w_v||U - Sv||^2 + lambda*trace(F'*Lu*F)}
% s.t Sv>=0, 1^T*Sv_i=1, U>=0, 1^T*Ui=1, F'*F=I
%
function [y, U, S0, S0_initial, F, evs] = GMC(X, c, lambda, normData)
%% input:
% X{}: multi-view dataset, each cell is a view, each column is a data point
% c: cluster number
% lambda: parameter (default 1)
%% output:
% S0: similarity-induced graph (SIG) matrix for each view
% y: the final clustering result, i.e., cluster indicator vector
% U: the learned unified matrix
% F: the embedding representation
% evs: eigenvalues of learned graph Laplacian matrix

NITER = 20;
zr = 10e-11;
pn = 15; % number of neighbours for constructS_PNG
islocal = 1; % only update the similarities of neighbors if islocal=1
if nargin < 3
    lambda = 1;
end;
if nargin < 4
    normData = 1;
end;

num = size(X{1},2); % number of instances
m = length(X); % number of views
%% Normalization: Z-score
if normData == 1
    for i = 1:m
        for  j = 1:num
            normItem = std(X{i}(:,j));
            if (0 == normItem)
                normItem = eps;
            end;
            X{i}(:,j) = (X{i}(:,j)-mean(X{i}(:,j)))/(normItem);
        end;
    end;
end;

%% initialize S0: Constructing the SIG matrices
S0 = cell(1,m);
for i = 1:m
    [S0{i}, ~] = InitializeSIGs(X{i}, pn, 0);
end;
S0_initial = S0;

%% initialize U, F and w
U = zeros(num);
for i = 1:m
    U = U + S0{i};
end;
U = U/m;
for j = 1:num
    U(j,:) = U(j,:)/sum(U(j,:));
end;
% % choose the top-k neighbors
% [~, ids] = sort(U,2,'descend');
% ts = zeros(num);
% for i =1:num
%     ts(i,ids(i,1:pn)) = U(i,ids(i,1:pn));
% end
% for j = 1:num
%     ts(j,:) = ts(j,:)/sum(ts(j,:));
% end
% U = ts;

sU = (U+U')/2;
D = diag(sum(sU));
L = D - sU;
[F, ~, evs]=eig1(L, c, 0);

w = ones(1,m)/m;

idxx = cell(1,m);
ed = cell(1,m);
for v = 1:m
    ed{v} = L2_distance_1(X{v}, X{v});
    [~, idxx{v}] = sort(ed{v}, 2); % sort each row
end;

%%  update ...
for iter = 1:NITER
    % update S^v
    for v = 1:m
        S0{v} = zeros(num);
        for i = 1:num
            id = idxx{v}(i,2:pn+2);
            di = ed{v}(i, id);
            numerator = di(pn+1)-di+2*w(v)*U(i,id(:))-2*w(v)*U(i,id(pn+1));
            denominator1 = pn*di(pn+1)-sum(di(1:pn));
            denominator2 = 2*w(v)*sum(U(i,id(1:pn)))-2*pn*w(v)*U(i,id(pn+1));
            S0{v}(i,id) = max(numerator/(denominator1+denominator2+eps),0);
        end;
%         for j = 1:num
%             normItem = sum(S0{v}(j,:));
%             if normItem == 0
%                 normItem = eps;
%             end;
%             S0{v}(j,:) = S0{v}(j,:)/normItem;
%         end;
    end;
    % update w
    for v = 1:m
        US = U - S0{v};
        distUS = norm(US, 'fro')^2;
        if distUS == 0
            distUS = eps;
        end;
        w(v) = 0.5/sqrt(distUS);
    end;
    % disp(['weights: ',num2str(w)]);
    % update U
    dist = L2_distance_1(F',F');
    U = zeros(num);
    for i=1:num
        idx = zeros();
        for v = 1:m
            s0 = S0{v}(i,:);
            idx = [idx,find(s0>0)];
        end;
        idxs = unique(idx(2:end));
        if islocal == 1
            idxs0 = idxs;
        else
            idxs0 = 1:num;
        end;
        for v = 1:m
            s1 = S0{v}(i,:);
            si = s1(idxs0);
            di = dist(i,idxs0);
            mw = m*w(v);
            lmw = lambda/mw;
            q(v,:) = si-0.5*lmw*di;
        end;
        U(i,idxs0) = SloutionToP19(q,m);
        clear q;
    end;
%         % choose the top-k neighbors
%         [~, ids] = sort(U,2,'descend');
%         ts = zeros(num);
%         for i =1:num
%             ts(i,ids(i,1:pn)) = U(i,ids(i,1:pn));
%         end
%         for j = 1:num
%             ts(j,:) = ts(j,:)/sum(ts(j,:));
%         end
%         sU = ts;
    % update F
    sU = U;
    sU = (sU+sU')/2;
    D = diag(sum(sU));
    L = D-sU;
    F_old = F;
    [F, ~, ev]=eig1(L, c, 0, 0);
    evs(:,iter+1) = ev;
    % update lambda and the stopping criterion
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > zr
        lambda = 2*lambda;
    elseif fn2 < zr
        lambda = lambda/2;
        F = F_old;
    else
        disp(['iter = ',num2str(iter),' lambda:',num2str(lambda)]);
        break;
    end;
end;
%% generating the clustering result
[clusternum, y]=graphconncomp(sparse(sU)); y = y';
if clusternum ~= c
    fprintf('Can not find the correct cluster number: %d\n', c)
end; 


