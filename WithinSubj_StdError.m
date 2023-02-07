function [SEM] = WithinSubj_StdError(Y_ij)
%factors on rows, (frequencies)
%subjects on columns (observations)
Y_i = mean(Y_ij,'omitnan');
Y_all = mean(Y_ij(:),'omitnan');
W_ij = zeros(size(Y_ij));
SEM_j = zeros(size(Y_ij,1),1);
W_j = SEM_j;
for i = 1 : size(Y_ij,1)
    for j = 1 : size(Y_ij,2)
        W_ij(i,j) = Y_ij(i,j) - (Y_i(i) - Y_all);        
    end
end
for i = 1 : size(Y_ij,1)
  SEM(i) = sqrt(mean((W_ij(i,:) -  mean(W_ij(i,:),'omitnan')).^2,'omitnan'))/ (size(W_ij,2) - 1);
end