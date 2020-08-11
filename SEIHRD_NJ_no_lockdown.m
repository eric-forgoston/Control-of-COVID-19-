%clear all;

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

par.rho=repmat(0.0397,1,151);

for j=1:17
for k=1:17
Crec(j,k)=0.5*(Contact.Cf(j,k)+Contact.Cf(k,j).*age.N(j)./age.N(k));
end
end
M=zeros(size(par.rho));
for l=1:151
R=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec;
M(1,l)=max(eig(R));
end
M;

ft=linspace(0,150,151);
[tf,yf]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[0:1:150],IC);

[p,q]=size(yf);

dlmwrite("herd.csv",yf);                 


Stotal=zeros(p,1);
Etotal=zeros(p,1);
Istotal=zeros(p,1);
Iatotal=zeros(p,1);
Htotal=zeros(p,1);
Rtotal=zeros(p,1);
Dtotal=zeros(p,1);

for k=1:17
    Stotal=Stotal+yf(:,0+k);
    Etotal=Etotal+yf(:,17+k);
    Istotal=Istotal+yf(:,34+k);
    Iatotal=Iatotal+yf(:,51+k);
    Htotal=Htotal+yf(:,68+k);
    Rtotal=Rtotal+yf(:,85+k);
    Dtotal=Dtotal+yf(:,102+k);
    
end

figure(1)
%subplot(3,1,1)
%plot(tf,Istotal,'k','LineWidth',3)
hold on
p1=plot(tf,yf(:,35),'Color',[.3686 .3098 .6353],'LineWidth',2)
p2=plot(tf,yf(:,36),'Color',[.2189 .4490 .7257],'LineWidth',2)
p3=plot(tf,yf(:,37),'Color',[.2170 .5938 .7268],'LineWidth',2)
p4=plot(tf,yf(:,38),'Color',[.3689 .7402 .6518],'LineWidth',2)
p5=plot(tf,yf(:,39),'Color',[.5331 .8193 .6450],'LineWidth',2)
p6=plot(tf,yf(:,40),'Color',[.7033 .8802 .6404],'LineWidth',2)
p7=plot(tf,yf(:,41),'Color',[.8663 .9507 .6030],'LineWidth',2)
p8=plot(tf,yf(:,42),'Color',[.9246 .9590 .6326],'LineWidth',2)
p9=plot(tf,yf(:,43),'Color',[.9500 .9500 .7115],'LineWidth',2)
p10=plot(tf,yf(:,44),'Color',[.9856 .9147 .6220],'LineWidth',2)
p11=plot(tf,yf(:,45),'Color',[.9957 .8434 .5025],'LineWidth',2)
p12=plot(tf,yf(:,46),'Color',[.9930 .7105 .3983],'LineWidth',2)
p13=plot(tf,yf(:,47),'Color',[.9804 .5539 .3044],'LineWidth',2)
p14=plot(tf,yf(:,48),'Color',[.9485 .4019 .2647],'LineWidth',2)
p15=plot(tf,yf(:,49),'Color',[.8735 .2911 .3024],'LineWidth',2)
p16=plot(tf,yf(:,50),'Color',[.7673 .1603 .3024],'LineWidth',2)
p17=plot(tf,yf(:,51),'Color',[.6196 .0039 .2588],'LineWidth',2)




ylabel('I_s','FontSize',16)
legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17],'[0-5]','[5-10]','[10-15]','[15-20]','[20-25]','[25-30]','[30-35]','[35-40]','[40-45]','[45-50]','[50-55]','[55-60]','[60-65]','[65-70]','[70-75]','[75-80]','[80+]','Location','NorthEast') 
xticks([0 50 100 150 200])
xticklabels({'2020-03-04','2020-04-23','2020-06-12','2020-08-01','2020-09-20'})
xtickangle(45)
box on

axes('Position',[.46 .6 .29 .29])
box on
plot(tf,Istotal,'k','LineWidth',3)
ylabel('I_s','FontSize',16)
xticks([0 50 100 150 200])
xticklabels({'2020-03-04','2020-04-23','2020-06-12','2020-08-01','2020-09-20'})
xtickangle(45)

figure(2)
hold on
p1=plot(tf,yf(:,103),'Color',[.3686 .3098 .6353],'LineWidth',2)
p2=plot(tf,yf(:,104),'Color',[.2189 .4490 .7257],'LineWidth',2)
p3=plot(tf,yf(:,105),'Color',[.2170 .5938 .7268],'LineWidth',2)
p4=plot(tf,yf(:,106),'Color',[.3689 .7402 .6518],'LineWidth',2)
p5=plot(tf,yf(:,107),'Color',[.5331 .8193 .6450],'LineWidth',2)
p6=plot(tf,yf(:,108),'Color',[.7033 .8802 .6404],'LineWidth',2)
p7=plot(tf,yf(:,109),'Color',[.8663 .9507 .6030],'LineWidth',2)
p8=plot(tf,yf(:,110),'Color',[.9246 .9590 .6326],'LineWidth',2)
p9=plot(tf,yf(:,111),'Color',[.9500 .9500 .7115],'LineWidth',2)
p10=plot(tf,yf(:,112),'Color',[.9856 .9147 .6220],'LineWidth',2)
p11=plot(tf,yf(:,113),'Color',[.9957 .8434 .5025],'LineWidth',2)
p12=plot(tf,yf(:,114),'Color',[.9930 .7105 .3983],'LineWidth',2)
p13=plot(tf,yf(:,115),'Color',[.9804 .5539 .3044],'LineWidth',2)
p14=plot(tf,yf(:,116),'Color',[.9485 .4019 .2647],'LineWidth',2)
p15=plot(tf,yf(:,117),'Color',[.8735 .2911 .3024],'LineWidth',2)
p16=plot(tf,yf(:,118),'Color',[.7673 .1603 .3024],'LineWidth',2)
p17=plot(tf,yf(:,119),'Color',[.6196 .0039 .2588],'LineWidth',2)



ylabel('D','FontSize',16)
legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17],'[0-5]','[5-10]','[10-15]','[15-20]','[20-25]','[25-30]','[30-35]','[35-40]','[40-45]','[45-50]','[50-55]','[55-60]','[60-65]','[65-70]','[70-75]','[75-80]','[80+]','Location','NorthEast') 
xticks([0 50 100 150 200])
xticklabels({'2020-03-04','2020-04-23','2020-06-12','2020-08-01','2020-09-20'})
xtickangle(45)
box on

axes('Position',[.2 .6 .29 .29])
box on
plot(tf,Dtotal,'k','LineWidth',3)
ylabel('D','FontSize',16)
xticks([0 50 100 150 200])
xticklabels({'2020-03-04','2020-04-23','2020-06-12','2020-08-01','2020-09-20'})
xtickangle(45)
