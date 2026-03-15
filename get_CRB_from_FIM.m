function [CRB_tau,CRB_nu,CRB_tau_sum,CRB_nu_sum] = get_CRB_from_FIM(FIM,invD,P)
    %CRB_matrix = invD*inv(invD*FIM*invD)*invD;
    CRB_matrix = inv(invD*FIM*invD);  %% Transformed CRB
    CRB_tau = diag(CRB_matrix(1:P,1:P)); 
    CRB_nu = diag(CRB_matrix(P+1:2*P,P+1:2*P));
    CRB_tau_sum = sum(CRB_tau);
    CRB_nu_sum = sum(CRB_nu);
end