
%% This function plot the number of antennas vs the sumrate and the no of anteannas vs the interefrence at the PU
%% This version will check the interefernce threshold level at the PU
%% Here we are using mAlso this does not use optimal powe allocation schemes as well
disp('======================================================');
disp('Sumrate for TWRN relay selection in a massive MIMO environment under a Primary Interferere');disp(' ');
Num = input('input # samples for estimation (1000s)  : ')*1000;
Nr=8;
eta1=1/16;
eta2=1/16;
% eta_k1=1*[1 1 1.01 0.99 1.02 0.98 1.01 1.01];
% eta_k2=1*[1 1 0.99 1.01 0.98 1.02 1.01 0.99];

eta_k1=1*[1 1 1.01 0.99 1.02 0.98 1.01 1.01];
eta_k2=1*[1 1 0.99 1.01 0.98 1.02 1.01 0.99];

%eta_PU=1/8*[1.001 1.002 1.003 0.999 0.998 1.002 0.998 1.002 0.997];

eta_PU=1/8*[1.01 1.02 1.03 0.99 0.98 1.02 0.98 1.02 0.97];
N=8;
Er=10^(7/10);
I_T=10^(15/10);
delta=.100;

sigmaK=1;
sigmaIk=1;
sigmaN=1;
sigmaIn=3;

N_vec=[10:5:30 40:10:100 120:20:200 220:40:300];

%N_vec=[10:40:100 100:100:300 400:300:2200];

K_vec=[1 2 8];
final_vec=zeros(3,size(N_vec,2));
final_int1=zeros(3,size(N_vec,2));
final_int2=zeros(3,size(N_vec,2));
outage_vec=zeros(3,size(N_vec,2));

for ghj=1:3
    K=K_vec(ghj);
   
    for i=1:size(N_vec,2)
        count=0;
        N1=N_vec(i);
        N2=N1;
        
        Pr=Er;
        Max_Sumrate=0;
        Inter1m=0;
        NumExt=Num;
        for p=1:NumExt
            Rl_vec = [];
            gamma_min_vec = [];
            Det_min_vec = [];
            F1  = (randn(N,N1) + 1j*randn(N,N1))*sqrt(1/2)*sqrt(eta1);
            F2  = (randn(N,N2) + 1j*randn(N,N2))*sqrt(1/2)*sqrt(eta2);
            SumRates=zeros(1,K);
            InfRates=zeros(1,K);
            Inter1m=Inter1m+P1*trace(F1'*F1)+P2*trace(F2'*F2);
            for ll=1:K                
                eta11=eta_k1(ll);
                eta22=eta_k2(ll);
                
                E1=I_T*eta22/((1+delta)*(eta1*eta22+eta2*eta11)*N);
                E2=I_T*eta11/((1+delta)*(eta1*eta22+eta2*eta11)*N);
                P1=E1/N1;
                P2=E2/N2;
                
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
                SumRates(1,ll)=2*sum(min(sum_rate1,sum_rate2));
                InfRates(1,ll)=Pr*trace(Gk'*Gk);
            end
            SumRates(InfRates>=I_T)=0;
            if P1*trace(F1'*F1)+P2*trace(F2'*F2)>=I_T
                SumRates=zeros(1,K);
            end
            [ret p]=max(SumRates);
            if(ret==0)
                NumExt=NumExt+1;
                count=count+1;
            end
            Max_Sumrate=Max_Sumrate+max(SumRates);
        end
        outage_vec(ghj,i)=count/Num;
        final_vec(ghj,i)=Max_Sumrate/Num;
        fprintf('N= %d Sumrate = %g\n',N1,  Max_Sumrate/Num);
    end
    
end
Pout1=(gammainc(I_T/(eta_PU(1)*Er),N*Nr,'upper'));
Pout2=(gammainc(I_T/(eta_PU(1)*Er),N*Nr,'upper'))*(gammainc(I_T/(eta_PU(2)*Er),N*Nr,'upper'));
Pout3=1;
for i=1:8
    Pout3=(gammainc(I_T/(eta_PU(i)*Er),N*Nr(1),'upper'))*Pout3;
end


analy1=(1-Pout1)*Nr*log(1+I_T*eta_k1(1)*eta_k2(1)/((1+delta)*(eta1*eta_k2(1)+eta2*eta_k1(1))*N*(sigmaK^2+sigmaIk^2)*Nr));
analy2=(1-Pout2)*Nr*log(1+I_T*min(eta_k1(1:2))*min(eta_k2(1:2))/((1+delta)*(eta1*min(eta_k2(1:2))+eta2*min(eta_k1(1:2)))*N*(sigmaK^2+sigmaIk^2)*Nr));
analy3=(1-Pout3)*Nr*log(1+I_T*min(eta_k1(1:8))*min(eta_k2(1:8))/((1+delta)*(eta1*min(eta_k2(1:8))+eta2*min(eta_k1(1:8)))*N*(sigmaK^2+sigmaIk^2)*Nr));




fig1=figure;

plot(N_vec, outage_vec(1,:), 'o--','color',[0 0 0],'LineWidth',2);
hold on;
plot(N_vec, outage_vec(2,:), 'o--','color',[0 1 0],'LineWidth',2);
hold on;
plot(N_vec, outage_vec(3,:), 'o--','color',[1 0 0],'LineWidth',2);
hold on;

plot(N_vec, Pout1*ones(size(final_vec(3,:))), '-','color',[0 0 0],'LineWidth',2);
hold on;
plot(N_vec, Pout2*ones(size(final_vec(3,:))), '-','color',[0 1 0],'LineWidth',2);
hold on;
plot(N_vec, Pout3*ones(size(final_vec(3,:))), '-','color',[1 0 0],'LineWidth',2);
hold on;
grid on;
legend('Case 1 (K=1)','Case 2 (K=2)','Case 3 (K=8)','Asymptotic case 1','Asymptotic case 2','Asymptotic case 3');
ylabel('Outage probability')


figure
plot(N_vec, final_vec(1,:), 'o--','color',[0 0 0],'LineWidth',2);
hold on;
plot(N_vec, final_vec(2,:), 'o--','color',[0 1 0],'LineWidth',2);
hold on;
plot(N_vec, final_vec(3,:), 'o--','color',[1 0 0],'LineWidth',2);
hold on;
plot(N_vec, analy1*ones(size(final_vec(3,:))), '-.','color',[0 0 0],'LineWidth',2);
hold on;
plot(N_vec, analy2*ones(size(final_vec(3,:))), '-','color',[0 1 0],'LineWidth',2);
hold on;
plot(N_vec, analy3*ones(size(final_vec(3,:))), '-','color',[1 0 0],'LineWidth',2);
hold on;
grid on;
legend('Case 1 (K=1)','Case 2 (K=2)','Case 3 (K=8)','Asymptotic case 1','Asymptotic case 2','Asymptotic case 3');
ylabel('Achievable sum rate (bps/Hz)')






