function [F_RF,F_BB]= OMP_precoding()
global Ns Ntf H AT 
F_RF = [];
[~,~,V]= svd(H);
F_opt = V(:,1,Ns);
F_res = F_opt;
for i = Ntf
    y = AT'*F_res;
    k = find(diag(y*y')==max(diag(y*y')));
    F_RF(:,i) = AT(:,k);
    F_BB=(F_RF'*F_RF)^(-1)*F_RF'*F_opt;
    F_res = (F_opt-F_RF*F_BB)/norm(F_opt-F_RF*F_BB,'fro');
end
F_BB=sqrt(Ns)*F_BB/norm(F_RF*F_BB,'fro');