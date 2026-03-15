% isotropic transmission varying SNR and number of targets
clc,clear all,close all
cd ..
run("parameter_setting.m");

SNR_list = 10:1:20;
pt_list = 10.^(SNR_list/10);

CRB_tau_isotropy_list = zeros(P,length(SNR_list));
CRB_nu_isotropy_list = zeros(P,length(SNR_list));
CRB_tau1_isotropy_list = zeros(P,length(SNR_list));
CRB_nu1_isotropy_list = zeros(P,length(SNR_list));

for p = 1:P
    D = diag([1/Tsym*ones(1,p),Tsym*ones(1,p),ones(1,2*p)]);
    invD = inv(D);
    [V,T,dotV,dotT] = channel_initialization(nu(1:p),tau(1:p),N,M,Tsym,Ts);
    Hq_overall = calculate_Hq(p,Q,g(1:p),V,dotV,T,dotT,noise_var);
    for pt_index = 1:length(pt_list)
        pt = pt_list(pt_index);
        FIM_isotropy = get_FIM_from_Hqq(ones(Q,1)*pt/Q,Hq_overall);
        [CRB_tau_tmp,CRB_nu_tmp,CRB_tau_isotropy_list(p,pt_index),CRB_nu_isotropy_list(p,pt_index)] = get_CRB_from_FIM(FIM_isotropy,invD,p);
        CRB_tau1_isotropy_list(p,pt_index) = CRB_tau_tmp(1);
        CRB_nu1_isotropy_list(p,pt_index) = CRB_nu_tmp(1);
    end

    if p == 1
        tau_analytical_list = 3./(2*abs(g(1))^2*pi^2*Q*(M^2-1)*(pt_list/noise_var));
        nu_analytical_list = 3./(2*abs(g(1))^2*pi^2*Q*(N^2-1)*(pt_list/noise_var));
    end
end

cd main_CRB

save('CRB_SNR_isotropy_data.mat','CRB_tau1_isotropy_list','CRB_nu1_isotropy_list', ...
    'CRB_tau_isotropy_list','CRB_nu_isotropy_list','tau_analytical_list','nu_analytical_list','SNR_list','-mat');