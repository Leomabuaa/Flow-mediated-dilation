% The program to calculate the expected FMD% for each patients. Copyright to Beihang University. 
%% Input the measured informations (measured by ultrasound and blood pressure monitor)
d_base=3e-3 % The inner arterial diameter at the end of the cuff inflated
D_base=4e-3 % The outer arterial diameter at the end of the cuff inflated
P_base=2660 % The brachial artery blood pressure at the end of the cuff inflated
U_fmd=[0.10;0.40] % The average blood flow velocity at the process of pre-inflation and the maximum blood flow velocity at the process of post-inflation
P_fmd=[7980;13300]% The average brachial artery blood pressure at the process of pre-inflation and the maximum brachial artery blood pressure at the process of post-inflation
%% Input basic parameters
Dl=3.3e-9 % Diffusivity in lumen 
Dw=8.48e-10 % Diffusivity in tissue
Kery=23 % Consumption in lumen
Kw=0.01 % Consumption in tissue
L=0.05 % Axial length in (m)
Z=0.2 % Axial position
miu=3.5e-3 % Viscousity of blood (newtonian)
niu=0.3 % Possion ratio
epi=1 % endothelial functionality
T=2e-6 % endothelial thickness
%% Derived parameters 
for i=1:1:2
    tao(i,1)=(4*miu*U_fmd(i,1))/d_base
    R(i,1)=2.13+(457.5*(tao(i,1)/(tao(i,1)+3.5)))
    E1=sqrt(1+(Kw*d_base^2)/Dw)
    E2=-1+E1
    E3=-1-E1
    E4=D_base/(2*d_base)
    E5=Dw/Dl
    E6(i,1)=(Kery*L)/U_fmd(i,1)
    E7(i,1)=exp(-E6(i,1)*Z)
    E8=d_base/Dl
    B25=(1+E1)*exp(E4*E3)*exp(0.5*E1)/(E2*exp(E4*E2))
    B26=exp(-0.5*E1)
    G5=E2*exp(E4*E2)
    G6=E3*exp(E4*E3)
    G9=-G6/G5
    L1(i,1)=(B25+B26)/E7(i,1)
    R1(i,1)= E8*R(i,1)*exp(0.5)*T*epi/E7(i,1)
    R2(i,1)=(E5*(1-(1/E5)+E1)*exp(-0.5*E1))/E7(i,1)
end
%% solve equations to get A and bw
for i=1:1:2
    syms x y
    [sol_x, sol_y] = vpasolve([((2*x+1)*0.000001/(4*x))*exp(x)==L1(i,1)*y, 
        ((2*x+1)*0.000001/(4*x))*exp(x) ==(R1(i,1)-R2(i,1)*y)/(2*x)],[x,y],[17;20000])
    sym_bw(i,1)=sol_y
end
for i=1:1:2
    bw(i,1)=double(sym_bw(i,1))
    aw(i,1)=G9*bw(i,1)
    c_mean(i,1)=((aw(i,1)*(exp(E4*E2)-exp(0.5*E2))/E2)+(bw(i,1)*(exp(E4*E3)-exp(0.5*E3))/E3))/(E4-0.5)
    E(i,1)=9.8e5*exp(-0.45*(c_mean(i,1)-1))+1e5;
end
%% Calculate the arterial deformation
dP=P_fmd-P_base
for i=1:1:2
   dd(i,1)=2*dP(i,1)*(1-niu*niu)*d_base^2*D_base/(E(i,1)*(D_base^2-d_base^2))
end
%% Calculate eFMD%
eFMD=(dd(2,1)-dd(1,1))*100/(dd(1,1)+d_base)
