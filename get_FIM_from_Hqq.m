function FIM = get_FIM_from_Hqq(rho_vec,Hq_overall)
    FIM = sum(Hq_overall .* reshape(rho_vec,1,1,[]), 3);
    % FIM = zeros(size(Hq_overall,[1,2]));
    % for qq = 1:length(rho_vec)
    %     FIM = FIM+rho_vec(qq)*Hq_overall(:,:,qq);
    % end
end