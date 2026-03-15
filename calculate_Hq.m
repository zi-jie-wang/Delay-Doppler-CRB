function Hq_overall = calculate_Hq(P,Q,g,V,dotV,T,dotT,noise_var)
    Hq_overall = calculate_FIM_slice(P,Q,g,V,dotV,T,dotT)*2*Q/noise_var;
end

function FIM = calculate_FIM_slice(P,Q,g,V,dotV,T,dotT)
    J_tautau = calculate_sub_block_slice(P,Q,g,g,V,V,dotT,dotT);
    J_taunu  = calculate_sub_block_slice(P,Q,g,g,V,dotV,T,dotT);
    J_nunu   = calculate_sub_block_slice(P,Q,g,g,dotV,dotV,T,T);
    J_taug   = calculate_sub_block_slice(P,Q,g,ones(1,P),V,V,T,dotT);
    J_nug    = calculate_sub_block_slice(P,Q,g,ones(1,P),dotV,V,T,T);
    J_gg     = calculate_sub_block_slice(P,Q,ones(1,P),ones(1,P),V,V,T,T);
    
    FIM = [
        real(J_tautau), real(J_taunu), real(J_taug), -imag(J_taug);    
        real(pagetranspose(J_taunu)), real(J_nunu), real(J_nug), -imag(J_nug);
        real(pagetranspose(J_taug)), real(pagetranspose(J_nug)), real(J_gg), -imag(J_gg);
        -imag(pagetranspose(J_taug)), -imag(pagetranspose(J_nug)), -imag(pagetranspose(J_gg)) real(J_gg);];
end

function J_tilde = calculate_sub_block_slice(P,Q,g1,g2,V1,V2,T1,T2)
    J_tilde = zeros(P,P,Q);
    for p = 1:P
        for pp = 1:P
            J_tilde(p,pp,:) = conj(g1(p))*g2(pp)* diag(kron(V1(:,:,p)'*V2(:,:,pp),T1(:,:,pp)*conj(T2(:,:,p))));
        end
    end
end