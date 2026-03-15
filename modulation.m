run('parameter_setting.m');

FN = dftmtx(N)/sqrt(N);
FM = dftmtx(M)/sqrt(M);
LAMBDA1 = diag(exp(-1j*2*pi*9/512*(0:1:M-1).^2));
LAMBDA2 = diag(exp(-1j*2*pi*(0)*(0:1:M-1).^2));
XI = zeros(M,M);
for m1 = 1:M
    for m2 = 1:M
        XI(m1,m2) = exp(1j*pi/M*(m1-m2)^2);
    end
end
XI = XI/sqrt(M)*exp(-1j*pi/4);

U = zeros(Q,Q,5);
U(:,:,1) = eye(Q); %%% OFDM
U(:,:,2) = kron(FN',FM); %%% OTFS;
U(:,:,3) = U(:,:,2); %%% ODDM
U(:,:,4) = kron(eye(N),FM*LAMBDA1*FM'*LAMBDA2);
U(:,:,5) = kron(eye(N),FM*XI);


waveform = {'OFDM','OTFS','ODDM','AFDM','OCDM'};