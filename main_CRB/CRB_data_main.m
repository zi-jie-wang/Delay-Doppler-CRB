clc,clear all,close all
cd ..
%run("parameter_setting.m");
run("modulation.m");

rng(2026);
data1 = kron((1:N)',ones(M,1))+kron(ones(N,1),(1:M)');
current_norm2 = sum(abs(data1).^2);
data1 = sqrt(pt*Q/current_norm2)*data1;

rho_vec = zeros(Q,size(U,3));
CRB_tau_list = zeros(size(U,3),P);
CRB_nu_list = zeros(size(U,3),P);

for u = 1:size(U,3)
    s = U(:,:,u)*data1;
    s = s * sqrt(pt*Q)/norm(s);
    rho_vec(:,u) = abs(s).^2/Q;
end

for p = 1:P
    D = diag([1/Tsym*ones(1,p),Tsym*ones(1,p),ones(1,2*p)]);
    invD = inv(D);
    [V,T,dotV,dotT] = channel_initialization(nu(1:p),tau(1:p),N,M,Tsym,Ts);
    Hq_overall = calculate_Hq(p,Q,g(1:p),V,dotV,T,dotT,noise_var);

    for u = 1:size(U,3)
        FIM_isotropy = get_FIM_from_Hqq(rho_vec(:,u),Hq_overall);
        [~,~,CRB_tau_list(u,p),CRB_nu_list(u,p)] = get_CRB_from_FIM(FIM_isotropy,invD,p);
    end
end

CRB_nu_list = sqrt(CRB_nu_list/(Tsym^2));
CRB_tau_list = sqrt(CRB_tau_list*(Tsym^2));

cd main_CRB;
save('CRB_data.mat','CRB_nu_list','CRB_tau_list','waveform','-mat');