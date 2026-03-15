function [V,T,dotV,dotT] = channel_initialization(nu,tau,N,M,Tsym,Ts)
    V = get_V(nu,N,Tsym);
    dotV = get_dotV(nu,N,Tsym);
    T = get_T(tau,M,Ts);
    dotT = get_dotT(tau,M,Ts);
end


function V = get_V(nu,N,Tsym)
    V = zeros(N,N,length(nu));
    for p = 1:length(nu)
        V(:,:,p) = diag(exp(1j*2*pi*nu(p)*(0:1:N-1)*Tsym));
    end
end

function T = get_T(tau,M,Ts)
    T = zeros(M,M,length(tau));
    for p = 1:length(tau)
        T(:,:,p) = diag(sinc((0:1:M-1)/M).^2.*exp(-1j*2*pi*tau(p)*(0:1:M-1)/(M*Ts)));
        %T(:,:,p) = diag(exp(-1j*2*pi*tau(p)*(0:1:M-1)/(M*Ts)));
    end
end

function dotV = get_dotV(nu,N,Tsym)
    dotV = zeros(N,N,length(nu));
    for p = 1:length(nu)
        dotV(:,:,p) = diag(1j*2*pi*(0:1:N-1)*Tsym.*exp(1j*2*pi*nu(p)*(0:1:N-1)*Tsym));
    end
end

function dotT = get_dotT(tau,M,Ts)
    dotT = zeros(M,M,length(tau));
    for p = 1:length(tau)
        dotT(:,:,p) = diag(-1j*(sinc((0:1:M-1)/M).^2)*2*pi.*(0:1:M-1)/(M*Ts).*exp(-1j*2*pi*tau(p)*(0:1:M-1)/(M*Ts)));
        %dotT(:,:,p) = diag(-1j*2*pi*(0:1:M-1)/(M*Ts).*exp(-1j*2*pi*tau(p)*(0:1:M-1)/(M*Ts)));
    end
end