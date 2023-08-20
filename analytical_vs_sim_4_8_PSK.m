
clc;
clear all;
close all;
for k=2:3;
    M=2^k;                       %size
    N = k*10^3;                    % number of symbols
    a = [0:M-1]*2*pi/M;            % phase values
    SNRdB  = [2:1:12];             % SNR range
    sdB  = SNRdB + 10*log10(k);
    
    % Mbinary to Gray code conversion
     b = [0:M-1];
     map =bitxor(b,floor(b/2));
     [tt ind] = sort(map);
    
    c = zeros(1,N);
    for i = 1:length(SNRdB)
        
        
        bits = rand(1,N*k,1)>0.5;                        % random 1's and 0's
        
        % binary to decimal
        bin2DecMatrix = ones(N,1)*(2.^((k-1):-1:0)) ;
        shape=  reshape(bits,k,N).';
        
        % decimal to binary
        G= (sum(shape.*bin2DecMatrix,2)).';
        
        % Gray code mapping
        dec = ind(G+1)-1; %
        ph= dec*2*pi/M;
        
        % modulation
        d= exp(1i*ph);
        s = d;
        
        % AWGN
        
        n = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
        
        % reciever
        r = s + 10^(-sdB(i)/20)*n;
        % demodulation
        
        e = angle(r);
        
        % phase
        
        e(e<0) = e(e<0) + 2*pi;
        c = 2*pi/M*round(e/(2*pi/M));
        
        c(c==2*pi) = 0;
        cd = round(c*M/(2*pi));
        
        % Decimal to Gray conversion
        
        f = map(cd+1);
        cb = dec2bin(f,k) ;
        
        cb = cb.';
        cb = cb(1:end).';
        cb = str2num(cb).' ;
        
        
        % errors
        Err(i) = size(find(bits- cb),2);
        
        
    end
    sBer=Err/(N*k);
    tBer =(1/k)*erfc(sqrt(k*10.^(SNRdB/10))*sin(pi/M));
    
    % plot
 
    semilogy(SNRdB,tBer,'rs-','LineWidth',2);
    hold on
    semilogy(SNRdB,sBer,'kx-','LineWidth',2);
    hold on
    
end
legend('QPSK Theoretical', 'QPSK Simulation','8PSK Theoretical','8PSK Simulation');
xlabel('SNR dB')
ylabel('Bit Error Rate')
title('BER VS SNR M-PSK')