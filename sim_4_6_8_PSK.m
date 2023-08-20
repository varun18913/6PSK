snr_db = 2:1:16;
A = 1;
num_sim = 10e5;

%Pe_simulation_QPSK
Pe_simulation = zeros(1,length(snr_db));
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
        elseif x>=0.75 && x
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

semilogy(snr_db, Pe_simulation,"g",Linewidth=2)
hold on






%Pe_simulation_6PSK
Pe_simulation_e = zeros(1,length(snr_db));
sxx = [A A/2 -A/2 -A -A/2 A/2];
syy = [0 sqrt(3)*A/2 sqrt(3)*A/2 0 -sqrt(3)*A/2 -sqrt(3)*A/2];
dss = zeros(1,6);

for j = 1:length(snr_db)
    snr_lin = 10^(0.1*snr_db(j));
    count = 0;
        
    for k = 1:1:num_sim
        n = sqrt(1/(2*snr_lin)).*[randn(1) randn(1)];
        x = rand;
        if x >= 0 && x < 0.16 
            sm = [sxx(1) syy(1)];
        elseif x >= 0.16 && x < 0.33
             sm = [sxx(2) syy(2)];
        elseif x >= 0.33 && x < 0.5
             sm = [sxx(3) syy(3)];
        elseif x>=0.5 && x < 0.66
             sm = [sxx(4) syy(4)];
        elseif x>=0.66 && x < 0.83
             sm = [sxx(5) syy(5)];
        elseif x>=0.83 && x < 1
             sm = [sxx(6) syy(6)];
        end
        
        r = sm + n;
    
        
        for i = 1:6 
            dss(i) = (r(1)-sxx(i)).^2 + (r(2)-syy(i)).^2;
        end
        
        if dss(1)==min(dss)
            shat = [sxx(1) syy(1)];
        elseif dss(2)==min(dss)
            shat = [sxx(2) syy(2)];
        elseif dss(3)==min(dss)
            shat = [sxx(3) syy(3)];
        elseif dss(4)==min(dss)
            shat = [sxx(4) syy(4)];
        elseif dss(5)==min(dss)
            shat = [sxx(5) syy(5)];
        elseif dss(6)==min(dss)
            shat = [sxx(6) syy(6)];
        end
        
        if(shat~=sm)
            count = count + 1;
        end 
        
    end
    
    Pe_simulation_e(j) = count/num_sim;
    
end

semilogy(snr_db, Pe_simulation_e,"b",Linewidth=2)
hold on

%Pe_simulation_8PSK
Pe_simulation_f = zeros(1,length(snr_db));
sxxx = [A A/sqrt(2) 0 -A/sqrt(2) -A -A/sqrt(2) 0 A/sqrt(2)];
syyy = [0 A/sqrt(2) A A/sqrt(2) 0 -A/sqrt(2) -A -A/sqrt(2)];
dsss = zeros(1,8); 

for j = 1:length(snr_db)
    snr_lin = 10^(0.1*snr_db(j));
    count = 0;
        
    for k = 1:1:num_sim
        n = sqrt(1/(2*snr_lin)).*[randn(1) randn(1)];
        x = rand;
        if x >= 0 && x < 0.125 
            sm = [sxxx(1) syyy(1)];
        elseif x >= 0.125 && x < 0.25
             sm = [sxxx(2) syyy(2)];
        elseif x >= 0.25 && x < 0.375
             sm = [sxxx(3) syyy(3)];
        elseif x>=0.375 && x < 0.5
             sm = [sxxx(4) syyy(4)];
        elseif x>=0.5 && x < 0.625
             sm = [sxxx(5) syyy(5)];
        elseif x>=0.625 && x < 0.75
             sm = [sxxx(6) syyy(6)];
        elseif x>=0.75 && x < 0.875
             sm = [sxxx(7) syyy(7)];
        elseif x>=0.875 && x < 1
             sm = [sxxx(8) syyy(8)];     
        end
    
        r = sm + n;
    
        
        for i = 1:8 
            dsss(i) = (r(1)-sxxx(i)).^2 + (r(2)-syyy(i)).^2;
        end
        
        if dsss(1)==min(dsss)
            shat = [sxxx(1) syyy(1)];
        elseif dsss(2)==min(dsss)
            shat = [sxxx(2) syyy(2)];
        elseif dsss(3)==min(dsss)
            shat = [sxxx(3) syyy(3)];
        elseif dsss(4)==min(dsss)
            shat = [sxxx(4) syyy(4)];
        elseif dsss(5)==min(dsss)
            shat = [sxxx(5) syyy(5)];
        elseif dsss(6)==min(dsss)
            shat = [sxxx(6) syyy(6)];
        elseif dsss(7)==min(dsss)
            shat = [sxxx(7) syyy(7)];
        elseif dsss(8)==min(dsss)
            shat = [sxxx(8) syyy(8)];  
        end
        
        if(shat~=sm)
            count = count + 1;
        end 
        
    end
    
    Pe_simulation_f(j) = count/num_sim;
    
    
end

semilogy(snr_db, Pe_simulation_f,"r",Linewidth=2)
hold on


% error_probability_derived=zeros(1,length(snr_db));
% for i=1:length(snr_db)
%     snr_linerar=10^(snr_db(i)/10);
%     error_probability_derived(i)=qfunc(sqrt(2*snr_linerar))+ ...
%         2*qfunc(sqrt(snr_linerar))*(1-qfunc(sqrt(2*snr_linerar)));
% end
% semilogy(snr_db,error_probability_derived,'*',Linewidth=1)
% hold on


xlabel('SNR [db]');
ylabel('Probability of Error');
title('Probability of Error vs SNR [db]');
legend('Pe:Simulation-QPSK', 'Pe:Simulation-6PSK','Pe:Simulation-8PSK','Pe_analytical-QPSK');



