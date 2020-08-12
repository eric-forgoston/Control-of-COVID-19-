clear all;

%setup
%age structure i=1..17 denoting age classes 0-5, 5-10, ..., 75-80, 80+.
%S(1) = y(1), S(2) = y(2), ..., S(17) = y(17)
%E(1) = y(18), E(2) = y(19), ..., E(16) = y(34)
%Is(1) = y(35), Is(2) = y(36), ..., Is(16) = y(51)
%Ia(1) = y(52), Ia(2) = y(53), ..., Ia(16) = y(68)
%H(1) = y(69), H(2) = y(70), ..., H(16) = y(85)
%R(1) = y(86), R(2) = y(87), ..., R(16) = y(102)
%D(1) = y(103), D(2) = y(104), ..., D(16) = y(119)

C=load('MUestimates_all_locations_2.csv');
Chome=load('MUestimates_home_2.csv');
Cwork=load('MUestimates_work_2.csv');
Cschool=load('MUestimates_school_2.csv');
Cother=load('MUestimates_other_locations_2.csv');

Cf=zeros(17,17);
Cf=C; Cf(1:16,17)=C(1:16,16); Cf(17,1:16)=C(16,1:16); Cf(17,17)=C(16,16);

Cfhome=zeros(17,17);
Cfhome=Chome; Cfhome(1:16,17)=Chome(1:16,16); Cfhome(17,1:16)=Chome(16,1:16); Cfhome(17,17)=Chome(16,16);

Cfwork=zeros(17,17);
Cfwork=Cwork; Cfwork(1:16,17)=Cwork(1:16,16); Cfwork(17,1:16)=Cwork(16,1:16); Cfwork(17,17)=Cwork(16,16);

Cfschool=zeros(17,17);
Cfschool=Cschool; Cfschool(1:16,17)=Cschool(1:16,16); Cfschool(17,1:16)=Cschool(16,1:16); Cfschool(17,17)=Cschool(16,16);

Cfother=zeros(17,17);
Cfother=Cother; Cfother(1:16,17)=Cother(1:16,16); Cfother(17,1:16)=Cother(16,1:16); Cfother(17,17)=Cother(16,16);

%Contact.Cf=Cf;
Contact.Cfhome=Cfhome;
Contact.Cfwork=Cfwork;
Contact.Cfschool=Cfschool;
Contact.Cfother=Cfother;

par.Ntot=8882190; 
age.N=[0.06*par.Ntot 0.06*par.Ntot 0.065*par.Ntot 0.065*par.Ntot 0.065*par.Ntot 0.065*par.Ntot 0.065*par.Ntot 0.065*par.Ntot 0.065*par.Ntot 0.065*par.Ntot 0.07*par.Ntot 0.07*par.Ntot 0.06*par.Ntot 0.06*par.Ntot 0.035*par.Ntot 0.035*par.Ntot 0.03*par.Ntot];

par.gamma=1/5;
age.f=[0.111 0.111 0.121 0.121 0.121 0.121 0.121 0.121 0.121 0.121 0.175 0.175 0.287 0.287 0.287 0.287 0.287];
par.sigmas=1/10;
par.sigmaa=1/4;
par.alphah=1/10.4; 
age.h=0.8*[0.182 0.055 0.055 0.055 0.068 0.068 0.139 0.139 0.139 0.139 0.251 0.251 0.251 0.512 0.512 0.512 0.617];
age.dh=2*[0.002 0 0 0 0.002 0.002 0.009 0.009 0.009 0.009 0.036 0.036 0.036 0.149 0.149 0.149 0.328];
age.d=2.5*[0.001 0 0 0 0.001 0.001 0.004 0.004 0.004 0.004 0.014 0.014 0.014 0.059 0.059 0.059 0.129];

n=17;
age.S=1:n;
age.E=(n+1):2*n;
age.Is=(2*n+1):3*n;
age.Ia=(3*n+1):4*n;
age.H=(4*n+1):5*n;
age.R=(5*n+1):6*n;
age.D=(6*n+1):7*n;

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Isic=25000; Eic=Isic/mean(age.f); Iaic=25000; 
Hic=0; Ric=0; Dic=0; Sic=par.Ntot-Eic-Isic-Iaic-Hic-Ric-Dic;
IC(age.S)=Sic*(age.N/par.Ntot);
IC(age.E)=Eic*(age.N/par.Ntot);
IC(age.Is)=Isic*(age.N/par.Ntot);
IC(age.Ia)=Iaic*(age.N/par.Ntot);
IC(age.H)=Hic*(age.N/par.Ntot);
IC(age.R)=Ric*(age.N/par.Ntot);
IC(age.D)=Dic*(age.N/par.Ntot);

Ctrlhome=eye(17,17);
Ctrlwork=eye(17,17);
Ctrlschool=eye(17,17);
Ctrlother=eye(17,17);
Contact.Cf=Ctrlhome*Contact.Cfhome+Ctrlwork*Contact.Cfwork+Ctrlschool*Contact.Cfschool+Ctrlother*Contact.Cfother;

par.rho=repmat(0.0397,1,13);

for j=1:17
for k=1:17
Crec(j,k)=0.5*(Contact.Cf(j,k)+Contact.Cf(k,j).*age.N(j)./age.N(k));
end
end
M=zeros(size(par.rho));
for l=1:13
R=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec;
M(1,l)=max(eig(R));
end
M;

ft=linspace(0,12,13);
[tf,yf]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[0:1:12],IC);

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p,q]=size(yf);
IC2(1,:)=yf(p,:);

par.rho=[0.0397:-0.0002:0.0387];

vw=0.1*ones(1,17); vo=0.1*ones(1,17); vh=1.0*ones(1,17);
Ctrlhome=diag(vh);
Ctrlwork=diag(vw);
Ctrlschool=zeros(17,17);
Ctrlother=diag(vo);
Contact.Cf=Ctrlhome*Contact.Cfhome+Ctrlwork*Contact.Cfwork+Ctrlschool*Contact.Cfschool+Ctrlother*Contact.Cfother;

for j=1:17
for k=1:17
Crec1(j,k)=0.5*(Contact.Cf(j,k)+Contact.Cf(k,j).*age.N(j)./age.N(k));
end
end
M2=zeros(size(par.rho));
for l=1:6
R1=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec1;
M2(1,l)=max(eig(R1));
end
M2;

ft=linspace(12,17,6);
[tf2,yf2]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[p-1:1:17],IC2);

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p2,q2]=size(yf2);
IC3(1,:)=yf2(p2,:);

par.rho=[0.0385:-0.0002:0.0301];

vw=0.1*ones(1,17); vo=0.1*ones(1,17); vh=1.0*ones(1,17);
Ctrlhome=diag(vh);
Ctrlwork=diag(vw);
Ctrlschool=zeros(17,17);
Ctrlother=diag(vo);
Contact.Cf=Ctrlhome*Contact.Cfhome+Ctrlwork*Contact.Cfwork+Ctrlschool*Contact.Cfschool+Ctrlother*Contact.Cfother;

for j=1:17
for k=1:17
Crec1(j,k)=0.5*(Contact.Cf(j,k)+Contact.Cf(k,j).*age.N(j)./age.N(k));
end
end
M3=zeros(size(par.rho));
for l=1:43
R1=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec1;
M3(1,l)=max(eig(R1));
end
M3;

ft=linspace(17,59,43);
[tf3,yf3]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1):1:59],IC3);

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p3,q3]=size(yf3);
IC4(1,:)=yf3(p3,:);

par.rho=repmat(0.0295,1,17);

vw=0.1*ones(1,17); vo=0.15*ones(1,17); vh=1.0*ones(1,17);
Ctrlhome=diag(vh);
Ctrlwork=diag(vw);
Ctrlschool=zeros(17,17);
Ctrlother=diag(vo);
Contact.Cf=Ctrlhome*Contact.Cfhome+Ctrlwork*Contact.Cfwork+Ctrlschool*Contact.Cfschool+Ctrlother*Contact.Cfother;

for j=1:17
for k=1:17
Crec2(j,k)=0.5*(Contact.Cf(j,k)+Contact.Cf(k,j).*age.N(j)./age.N(k));
end
end
M4=zeros(size(par.rho));
for l=1:17
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M4(1,l)=max(eig(R2));
end
M4;

ft=linspace(59,75,17);
[tf4,yf4]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1):1:75],IC4);

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p4,q4]=size(yf4);
IC5(1,:)=yf4(p4,:);

par.rho=repmat(0.0295,1,5);

vw=0.15*ones(1,17); vo=0.15*ones(1,17); vh=1.0*ones(1,17);
Ctrlhome=diag(vh);
Ctrlwork=diag(vw);
Ctrlschool=zeros(17,17);
Ctrlother=diag(vo);
Contact.Cf=Ctrlhome*Contact.Cfhome+Ctrlwork*Contact.Cfwork+Ctrlschool*Contact.Cfschool+Ctrlother*Contact.Cfother;

for j=1:17
for k=1:17
Crec2(j,k)=0.5*(Contact.Cf(j,k)+Contact.Cf(k,j).*age.N(j)./age.N(k));
end
end
M5=zeros(size(par.rho));
for l=1:5
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M5(1,l)=max(eig(R2));
end
M5;

ft=linspace(75,79,5);
[tf5,yf5]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1):1:79],IC5);

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p5,q5]=size(yf5);
IC6(1,:)=yf5(p5,:);

Yf=yf5;
YFIs=Yf(:,35:51);
sz=size(YFIs);
Is=zeros(sz(1),1);
for i=1:sz
Is(i)=sum(YFIs(i,:));
end
Isbase=max(Is);

par.rho=repmat(0.0295,1,5);

vw=0.15*ones(1,17); vo=0.2*ones(1,17); vh=1.0*ones(1,17);
Ctrlhome=diag(vh);
Ctrlwork=diag(vw);
Ctrlschool=zeros(17,17);
Ctrlother=diag(vo);
Contact.Cf=Ctrlhome*Contact.Cfhome+Ctrlwork*Contact.Cfwork+Ctrlschool*Contact.Cfschool+Ctrlother*Contact.Cfother;

for j=1:17
for k=1:17
Crec2(j,k)=0.5*(Contact.Cf(j,k)+Contact.Cf(k,j).*age.N(j)./age.N(k));
end
end
M6=zeros(size(par.rho));
for l=1:5
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M6(1,l)=max(eig(R2));
end
M6;

ft=linspace(79,83,5);
[tf6,yf6]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1):1:83],IC6);

%----------------------------------------------------------------------------------------------------------------------------------------------


[p6,q6]=size(yf6);
IC7(1,:)=yf6(p6,:);


Yf=yf6;
YFIs=Yf(:,35:51);
sz=size(YFIs);
Is=zeros(sz(1),1);
for i=1:sz
Is(i)=sum(YFIs(i,:));
end
Isratio=max(Is)/Isbase;
if Isratio < 0.1
lock=0;
else 
lock=1;
end

par.rho=repmat(0.0295,1,7);


vw=0.2*ones(1,17); vo=0.2*ones(1,17); vh=1.0*ones(1,17);
Ctrlhome=diag(vh);
Ctrlwork=diag(vw);
Ctrlschool=zeros(17,17);
Ctrlother=diag(vo);
Contact.Cf=Ctrlhome*Contact.Cfhome+Ctrlwork*Contact.Cfwork+Ctrlschool*Contact.Cfschool+Ctrlother*Contact.Cfother;

for j=1:17
for k=1:17
Crec2(j,k)=0.5*(Contact.Cf(j,k)+Contact.Cf(k,j).*age.N(j)./age.N(k));
end
end
M7=zeros(size(par.rho));
for l=1:7
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M7(1,l)=max(eig(R2));
end
M7;


ft=linspace(83,89,7);
[tf7,yf7]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1):1:89],IC7);


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

[p7,q7]=size(yf7);
stincr=(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1);

baseYF=[yf; yf2; yf3; yf4; yf5; yf6; yf7];
baseRo=[M M2 M3 M4 M5 M6 M7];
cycYF=baseYF;
cycRo=baseRo;

ICS(1,:)=yf7(p7,:);

vwi=0.2;
vsi=0;
voi=0.2;

s=89;
t=s+10;

for n1=1:190



par.rho=repmat(0.0295,1,10);


if lock 
vwi=vwi-0.1;
voi=voi-0.1;
vsi=0;
else 

if vwi < 1
vwi=vwi+0.1;
end

if voi < 1
voi=voi+0.1;
end

if s>181
vsi=1;
end

end

vw=vwi*ones(1,17); vo=voi*ones(1,17); vh=1.0*ones(1,17);vs=vsi*ones(1,17);
Ctrlhome=diag(vh);
Ctrlwork=diag(vw);
Ctrlschool=diag(vs);
Ctrlother=diag(vo);
Contact.Cf=Ctrlhome*Contact.Cfhome+Ctrlwork*Contact.Cfwork+Ctrlschool*Contact.Cfschool+Ctrlother*Contact.Cfother;

for j=1:17
for k=1:17
Crec2(j,k)=0.5*(Contact.Cf(j,k)+Contact.Cf(k,j).*age.N(j)./age.N(k));
end
end
Mincr=zeros(size(par.rho));
for l=1:10
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
Mincr(1,l)=max(eig(R2));
end
Mincr;

ft=linspace(s,t,10);
[tfincr,yfincr]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[stincr:1:t],ICS);

[pincr,qincr]=size(yfincr);
stincr=stincr+(pincr-1);
ICS(1,:)=yfincr(pincr,:);
cycYF=[cycYF; yfincr];
cycRo=[cycRo Mincr];

YFIs=yfincr(:,35:51);
sz=size(YFIs);
Is=zeros(sz(1),1);
for i=1:sz
Is(i)=sum(YFIs(i,:));
end
Isratio=max(Is)/Isbase;
if Isratio < 0.1
lock=0;
else
lock=1;
end

s=s+10;
t=s+10;

end

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dlmwrite("lockdowncyclical.csv",cycYF);                 
dlmwrite("lockdowncyclicalRo.csv",cycRo);

[sig,tau]=size(cycYF);


Iatotal=zeros(sig,1);
Istotal=zeros(sig,1);
Htotal=zeros(sig,1);
Dtotal=zeros(sig,1);


for k=1:17
    Iatotal=Iatotal+cycYF(:,51+k);
    Istotal=Istotal+cycYF(:,34+k);
    Htotal=Htotal+cycYF(:,68+k);
    Dtotal=Dtotal+cycYF(:,102+k);   
end

tfin=[0:1:2185]';

figure(1)


subplot(4,1,1)
plot(tfin,Iatotal,'k','LineWidth',3)
ylabel('I_a','FontSize',16)
xticks([10 365 730 1095 1460 1825 2190])
xticklabels({'03/14/20','03/14/21','03/14/22','03/14/23','03/14/24','03/14/25','03/14/26'})
set(gca,'Xticklabel',[])
xtickangle(45)
axis([10 2200 0 5000])

subplot(4,1,2)
plot(tfin,Istotal,'k','LineWidth',3)
ylabel('I_s','FontSize',16)
xticks([10 365 730 1095 1460 1825 2190])
xticklabels({'03/14/20','03/14/21','03/14/22','03/14/23','03/14/24','03/14/25','03/14/26'})
set(gca,'Xticklabel',[])
xtickangle(45)
axis([10 2200 0 2000])


subplot(4,1,3)
plot(tfin,Htotal,'k','LineWidth',3)
ylabel('H','FontSize',16)
xticks([10 365 730 1095 1460 1825 2190])
xticklabels({'03/14/20','03/14/21','03/14/22','03/14/23','03/14/24','03/14/25','03/14/26'})
set(gca,'Xticklabel',[])
xtickangle(45)
axis([10 2200 0 1000])


subplot(4,1,4)
plot(tfin,Dtotal,'k','LineWidth',3)
ylabel('D','FontSize',16)
xticks([10 365 730 1095 1460 1825 2190])
xticklabels({'2020-03-14','2021-03-14','2022-03-14','2023-03-14','2024-03-14','2025-03-14','2026-03-14'})
xtickangle(45)
axis([0 2200 0 20000])





