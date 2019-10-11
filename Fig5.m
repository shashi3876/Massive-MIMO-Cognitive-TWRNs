%% This function plot the number of antennas vs the sumrate and the no of anteannas vs the interefrence at the PU
%% This version will check the interefernce threshold level at the PU
%% Here we are using mAlso this does not use optimal powe allocation schemes as well
disp('======================================================');
disp('Sumrate for TWRN relay selection in a massive MIMO environment under a Primary Interferere');disp(' ');
Num = input('input # samples for estimation (1000s)  : ')*1000;
Nr_vec=[8 8 8 8 8 8 8 8];
eta1=1/16;
eta2=1/16;
eta_k1=2*[1 1 1 1 1 1 1 1];
eta_k2=2*[1 1 1 1 1 1 1 1];
eta_PU=1*[1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8];
N=6;
Er=10^(25/10);
I_T=10^(10/10);
delta=.100;
E1=I_T/((1+delta)*(eta1+eta2)*N)
E2=E1;
% ===========================================================================
% Simulation of BER
sigmaK=1;
sigmaIk=1;
sigmaN=1;
sigmaIn=1;
N_vec=[10:10:100 110:10:200 220:20:400];
K_vec=[1 4 8];
final_vec=zeros(3,size(N_vec,2));
for ghj=1:3
    K=K_vec(ghj);
    Sum_rate=[];
    for i=1:size(N_vec,2)
        N1=N_vec(i);
        N2=N1;
        P1=E1/N1;
        P2=E2/N2;
        Pr=Er/N1;        
        Max_Sumrate=0;
        Inter1m=0;      
        for p=1:Num
            F1  = (randn(N,N1) + 1j*randn(N,N1))*sqrt(1/2)*sqrt(eta1);
            F2  = (randn(N,N2) + 1j*randn(N,N2))*sqrt(1/2)*sqrt(eta2);
            SumRates=zeros(1,K);
            InfRates=zeros(1,K);            
            Inter1m=Inter1m+P1*trace(F1'*F1)+P2*trace(F2'*F2);
            for ll=1:K
                Nr=Nr_vec(ll);
                eta11=eta_k1(ll);
                eta22=eta_k2(ll);
                H1  = (randn(Nr,N1) + 1j*randn(Nr,N1))*sqrt(1/2)*sqrt(eta11);
                H2  = (randn(Nr,N2) + 1j*randn(Nr,N2))*sqrt(1/2)*sqrt(eta22);               
                Gk  = (randn(N,Nr) + 1j*randn(N,Nr))*sqrt(1/2)*sqrt(eta_PU(ll));
                              BN1=inv(H1*H1');
                BN2=inv(H2*H2');
                m1=sqrt(P1/trace(BN1));
                m2=sqrt(P2/trace(BN2));                
                Mk=sqrt(Pr/(m1^2+m2^2+sigmaK^2+sigmaIk^2));                
                det1=diag(BN1);
                det2=diag(BN2);
                gamma1=Mk^2*m2^2*ones(size(det1))./(Mk^2*(sigmaK^2+sigmaIk^2)+(sigmaIn^2+sigmaN^2)*det1);
                gamma2=Mk^2*m1^2*ones(size(det2))./(Mk^2*(sigmaK^2+sigmaIk^2)+(sigmaIn^2+sigmaN^2)*det2);
                sum_rate1=0.5*log(1+gamma1);
                sum_rate2=0.5*log(1+gamma2);
                SumRates(1,ll)=2*min(sum(sum(sum_rate1)),sum(sum(sum_rate2)));
                InfRates(1,ll)=Pr*trace(Gk'*Gk);
                %pause
            end
            SumRates(InfRates>I_T)=0;            
            if P1*trace(F1'*F1)+P2*trace(F2'*F2)>=I_T
                SumRates=zeros(1,K);
            end
            %==========================================================================
            [ret p]=max(SumRates);
            Max_Sumrate=Max_Sumrate+max(SumRates);  
        end        
        Sum_rate=[Sum_rate Max_Sumrate/Num];
        fprintf('N= %d Sumrate = %g\n',N1,  Max_Sumrate/Num);
    end
      final_vec(ghj,:)=Sum_rate;
end

analy3=1*max(Nr_vec)*log(1+Er*E2*eta11*eta22/((sigmaIn^2+sigmaN^2)*(E1*eta11+E2*eta22)+(sigmaK^2+sigmaIk^2)*max(Nr_vec)*(Er*eta11+sigmaIn^2+sigmaN^2)));


plot(N_vec, final_vec(1,:), 'o--','color',[0 0 0],'LineWidth',2);
hold on;
plot(N_vec, final_vec(2,:), 'o--','color',[0 1 0],'LineWidth',2);
hold on;
plot(N_vec, final_vec(3,:), 'o--','color',[1 0 0],'LineWidth',2);
hold on;
plot(N_vec, analy3*ones(size(final_vec(3,:))), 'x-','color',[0 0 0],'LineWidth',2);
hold on;
grid on;
legend('Case 1 (K=1)','Case 2 (K=4)','Case 3 (K=8)','Asymptotic value');
ylabel('Achievable sum rate (bps/Hz)')
xlabel('Number of antennas (N_1)')






