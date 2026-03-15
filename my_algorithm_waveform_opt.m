function [rho_vec,FIM_opt] = my_algorithm_waveform_opt(pt,P,Q,Hq_overall,invD,AAalpha,rho_vec_lastround,eta,gradiant_threshold_converge,function_threshold_converge)
    if rho_vec_lastround == 0
        rho_vec = pt/Q*ones(Q,1);
    else
        %%%% initialize rho_vec using last optimal rho_vec for adjacent alpha
        rho_vec = rho_vec_lastround;   
    end
    k = 0;
    tmp_grad = 1e-6;
    tmp_CRB_tau = 1;
    tmp_CRB_nu = 1;
    eye_rep = repmat(eye(size(Hq_overall,1)), [1 1 Q]);
    while true
        Jk = get_FIM_from_Hqq(rho_vec,Hq_overall);
        invJk = invD/(invD*Jk*invD)*invD;
        term1 = invJk*AAalpha*invJk;
        diag_only = pagemtimes(term1,Hq_overall) .* eye_rep;
        Deltak = squeeze(sum(diag_only, [1 2]));
        rho_vec = max(zeros(Q,1),rho_vec+eta*(Deltak));  
        rho_vec = rho_vec*pt/sum(rho_vec);
        FIM_tmp = get_FIM_from_Hqq(rho_vec,Hq_overall);%%%
        [~,~,CRB_tau,CRB_nu] = get_CRB_from_FIM(FIM_tmp,invD,P); %%%
        k = k+1;
        if abs(vecnorm(Deltak)/tmp_grad - 1) < gradiant_threshold_converge && abs(CRB_tau/tmp_CRB_tau-1) < function_threshold_converge && abs(CRB_nu/tmp_CRB_nu-1) < function_threshold_converge
            break
        end
        tmp_grad = vecnorm(Deltak);
        tmp_CRB_nu = CRB_nu;
        tmp_CRB_tau = CRB_tau;
    end
    FIM_opt = FIM_tmp;%get_FIM_from_Hqq(rho_vec,Hq_overall);
end