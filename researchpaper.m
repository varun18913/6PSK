%QPSK with differential encoding
snr_db = 2:1:12;
Pe_simulationx = zeros(1,length(snr_db));
for i = 1:length(snr_db)
    Pe_simulationx(i) = cal(snr_db(i));
end

semilogy(snr_db,Pe_simulationx,'-or',Linewidth=2)
xlabel('SNR [db]');
ylabel('Probability of Error');
title('Probability of Error vs SNR [db]');

hold on

%QPSK w/o differential encoding

A = 1;
num_sim = 10e5;
Pe_simulation = zeros(1,length(snr_db));
xlabel('SNR [db]');
ylabel('Probability of Error');
title('Probability of Error vs SNR [dB]');
sx = [A 0 -A 0];
sy = [0 A 0 -A];
ds = zeros(1,4);
for j = 1:length(snr_db)
    snr_lin = 10^(0.1*snr_db(j));
    count = 0;
        
    for k = 1:1:num_sim
        n = sqrt(1/(2*snr_lin)).*[randn(1) randn(1)];
        x = rand;
        if x >= 0 && x < 0.25 
            sm = [sx(1) sy(1)];
        elseif x >= 0.25 && x < 0.5
             sm = [sx(2) sy(2)];
        elseif x >= 0.5 && x < 0.75
             sm = [sx(3) sy(3)];
        elseif x>=0.75 && x <=1
             sm = [sx(4) sy(4)];
        end
    
        r = sm + n;
    
        
        for i = 1:4 
            ds(i) = (r(1)-sx(i)).^2 + (r(2)-sy(i)).^2;
        end
        
        if ds(1)==min(ds)
            shat = [sx(1) sy(1)];
        elseif ds(2)==min(ds)
            shat = [sx(2) sy(2)];
        elseif ds(3)==min(ds)
            shat = [sx(3) sy(3)];
        elseif ds(4)==min(ds)
            shat = [sx(4) sy(4)];
        end
        
        if(shat~=sm)
            count = count + 1;
        end 
        
    end
    
    Pe_simulation(j) = count/num_sim;
    
    
end

semilogy(snr_db, Pe_simulation,'b',Linewidth=2)
legend('Pe:simulation(diff)','Pe:simulation(w/o diff)')


%differential encoding
function y = dpske(A,M,mapping,sequence)
    k = log2(M);
    N = length(sequence);

    remainder = rem(N,k);
    if (remainder~=0)
        for i = N+1:N+k-remainder
            sequence(i) = 0;
        end
        N = N+k-remainder;
    end
    theta = 0;
    for i = 1:k:N
        index = 0;
        for j = i:i+k-1
            index = 2*index + sequence(j);
        end
        index = index + 1;
        theta = mod(2*pi*mapping(index)/M + theta, 2*pi);
        y((i+k-1)/k,1)=sqrt(A)*cos(theta);
        y((i+k-1)/k,2)=sqrt(A)*sin(theta);
    end

end




%differential encoding BER for M = 4
function y = cal(snr_e)
    N = 10000;
    A = 1;
    snr = 10^(snr_e/10);
    for i = 1:2*N
        temp = rand;
        if(temp<0.5)
            dsource(i) = 0;
        else
            dsource(i) = 1;
        end
    end
    
    mapping = [0 1 3 2];
    M = 4;
    diff_enc_output = dpske(A,M,mapping,dsource);
    
    for i = 1:N
        n = sqrt(1/(2*snr))*randn(1);
        r(i,:) = diff_enc_output(i,:)+n;
    end
    
    numoferr=0;
    prev_theta = 0;
    for i = 1:N
        theta = angle(r(i,1)+1i*r(i,2));
        delta_theta = mod(theta-prev_theta,2*pi);
        if((delta_theta<pi/4) || delta_theta>7*pi/4)
            decis = [ 0 0];
        elseif(delta_theta<3*pi/4)
            decis = [ 0 1];
        elseif(delta_theta<5*pi/4)
            decis = [ 1 1];
        elseif(delta_theta<7*pi/4)
            decis = [ 1 0 ];
        end
       
        prev_theta = theta;
        
        if((decis(1) ~= dsource(2*i-1) || decis(2) ~= dsource(2*i)))
            numoferr = numoferr+1;
        end
    end
    y = numoferr/N;
end





