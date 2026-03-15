clc,clear all,close all
cd ..
run("parameter_setting.m");

[V,T,dotV,dotT] = channel_initialization(nu,tau,N,M,Tsym,Ts);
Hq_overall = calculate_Hq(P,Q,g,V,dotV,T,dotT,noise_var);
save('Hq_default.mat','Hq_overall','-mat');

D = diag([1/Tsym*ones(1,P),Tsym*ones(1,P),ones(1,2*P)]);
invD = inv(D);

alpha_list = linspace(0,1,21);

CRB_tau_list = zeros(1,length(alpha_list));
CRB_nu_list = zeros(1,length(alpha_list));
CRB_tau_list_CVX = zeros(1,length(alpha_list));
CRB_nu_list_CVX = zeros(1,length(alpha_list));
rho_vec = zeros(Q,length(alpha_list));
rho_vec_last_round = 0;

for a = length(alpha_list):-1:1
    a
    alpha = alpha_list(a);
    AAalpha = diag([alpha/Tsym^2*ones(1,P),(1-alpha)*Tsym^2*ones(1,P),zeros(1,2*P)]);
    [rho_vec_last_round,FIM_opt] = my_algorithm_waveform_opt(pt,P,Q,Hq_overall,invD,AAalpha,rho_vec_last_round,300,1e-7,0.005);
    rho_vec(:,a) = rho_vec_last_round;
    [~,~,CRB_tau_list(a),CRB_nu_list(a)] = get_CRB_from_FIM(FIM_opt,invD,P);
    %[CRB_tau_list_CVX(a),CRB_nu_list_CVX(a)] = waveform_opt_CVX(pt,P,Q,invD,AAalpha,V,T,dotV,dotT);
end

plot(CRB_tau_list,CRB_nu_list);
hold on
FIM_isotropy = get_FIM_from_Hqq(ones(Q,1)*pt/Q,Hq_overall);
[~,~,CRB_tau_isotropy,CRB_nu_isotropy] = get_CRB_from_FIM(FIM_isotropy,invD,P);
scatter(CRB_tau_isotropy,CRB_nu_isotropy);

%%%% single-target ONLY!!!
% if P == 1
%     tt = 3/(2*abs(g(1))^2*pi^2*Q*(M^2-1)*(pt/noise_var));
%     vv = 3/(2*abs(g(1))^2*pi^2*Q*(N^2-1)*(pt/noise_var));
%     scatter(tt,vv)
% end

xlabel('$\varepsilon(\mathbf{\tau})$','Interpreter','latex');
cd Pareto_Front;
save('Data1.mat','alpha_list','rho_vec','CRB_tau_isotropy','CRB_nu_isotropy','CRB_tau_list','CRB_nu_list','-mat')