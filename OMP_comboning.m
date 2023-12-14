function [W_RF.W_BB]= OMP_combining(F_RF,F_BB)
global Ns Nr Nrf H AR Vn
W_RF = [];
W_MMSE= ((F_BB'*F_RF'*H'*H*F_RF*F_BB+Vn*Ns*eye(Ns))^(-1)*F_BB'*F_RF'*H')';
W_res = W_opt;
n = 1/Ns*H*F_RF*F_BB*F_BB'*F_RF'*H'+Vn*eye(Nr);
for i = Nrf
    y = AR'*n*W_res;
    k = find(diag(y*y')==max(diag(y*y')));
    W_RF(:,i) = AR(:,k);
    W_BB=(W_RF'*n*W_RF)^(-1)*W_RF'*n*W_MMSE;
    W_res = (W_MMSE-W_RF*W_BB)/norm(W_MMSE-W_RF*W_BB,'fro');
end
