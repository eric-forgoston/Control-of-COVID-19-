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

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p6,q6]=size(yf6);
IC7(1,:)=yf6(p6,:);

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


[p7,q7]=size(yf7);
IC8(1,:)=yf7(p7,:);


par.rho=repmat(0.0295,1,92);


vw=1.0*ones(1,17); vo=1.0*ones(1,17); vh=1.0*ones(1,17);vs=0.0*ones(1,17);
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
M8=zeros(size(par.rho));
for l=1:92
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M8(1,l)=max(eig(R2));
end
M8;


ft=linspace(89,180,92);
[tf8,yf8]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1):1:180],IC8);


[p8,q8]=size(yf8);
IC9(1,:)=yf8(p8,:);


par.rho=repmat(0.0295,1,360);


vw=1.0*ones(1,17); vo=1.0*ones(1,17); vh=1.0*ones(1,17);vs=1.0*ones(1,17);
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
M9=zeros(size(par.rho));
for l=1:360
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M9(1,l)=max(eig(R2));
end
M9;

ft=linspace(180,539,360);
[tf9,yf9]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1)+(p8-1):1:539],IC9);

[p9,q9]=size(yf9);
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

YF=[yf; yf2; yf3; yf4; yf5; yf6; yf7; yf8];
Ro=[M M2 M3 M4 M5 M6 M7 M8];        
dlmwrite("lockdownJun1.csv",YF);                 
dlmwrite("lockdownRoJun1.csv",Ro);

Stotal=zeros(p,1);
Stotal2=zeros(p2,1);
Stotal3=zeros(p3,1);
Stotal4=zeros(p4,1);
Stotal5=zeros(p5,1);
Stotal6=zeros(p6,1);
Stotal7=zeros(p7,1);
Stotal8=zeros(p8,1);
Stotal9=zeros(p9,1);
Etotal=zeros(p,1);
Etotal2=zeros(p2,1);
Etotal3=zeros(p3,1);
Etotal4=zeros(p4,1);
Etotal5=zeros(p5,1);
Etotal6=zeros(p6,1);
Etotal7=zeros(p7,1);
Etotal8=zeros(p8,1);
Etotal9=zeros(p9,1);
Istotal=zeros(p,1);
Istotal2=zeros(p2,1);
Istotal3=zeros(p3,1);
Istotal4=zeros(p4,1);
Istotal5=zeros(p5,1);
Istotal5=zeros(p5,1);
Istotal6=zeros(p6,1);
Istotal7=zeros(p7,1);
Istotal8=zeros(p8,1);
Istotal9=zeros(p9,1);
Iatotal=zeros(p,1);
Iatotal2=zeros(p2,1);
Iatotal3=zeros(p3,1);
Iatotal4=zeros(p4,1);
Iatotal5=zeros(p5,1);
Iatotal6=zeros(p6,1);
Iatotal7=zeros(p7,1);
Iatotal8=zeros(p8,1);
Iatotal9=zeros(p9,1);
Htotal=zeros(p,1);
Htotal2=zeros(p2,1);
Htotal3=zeros(p3,1);
Htotal4=zeros(p4,1);
Htotal5=zeros(p5,1);
Htotal5=zeros(p5,1);
Htotal6=zeros(p6,1);
Htotal7=zeros(p7,1);
Htotal8=zeros(p8,1);
Htotal9=zeros(p9,1);
Rtotal=zeros(p,1);
Rtotal2=zeros(p2,1);
Rtotal3=zeros(p3,1);
Rtotal4=zeros(p4,1);
Rtotal5=zeros(p5,1);
Rtotal6=zeros(p6,1);
Rtotal7=zeros(p7,1);
Rtotal8=zeros(p8,1);
Rtotal9=zeros(p9,1);
Dtotal=zeros(p,1);
Dtotal2=zeros(p2,1);
Dtotal3=zeros(p3,1);
Dtotal4=zeros(p4,1);
Dtotal5=zeros(p5,1);
Dtotal5=zeros(p5,1);
Dtotal6=zeros(p6,1);
Dtotal7=zeros(p7,1);
Dtotal8=zeros(p8,1);
Dtotal9=zeros(p9,1);
for k=1:17
    Stotal=Stotal+yf(:,0+k);
    Stotal2=Stotal2+yf2(:,0+k);
    Stotal3=Stotal3+yf3(:,0+k);
    Stotal4=Stotal4+yf4(:,0+k);
    Stotal5=Stotal5+yf5(:,0+k);
    Stotal6=Stotal6+yf6(:,0+k);
    Stotal7=Stotal7+yf7(:,0+k);
    Stotal8=Stotal8+yf8(:,0+k);
    Stotal9=Stotal9+yf9(:,0+k);
    Etotal=Etotal+yf(:,17+k);
    Etotal2=Etotal2+yf2(:,17+k);
    Etotal3=Etotal3+yf3(:,17+k);
    Etotal4=Etotal4+yf4(:,17+k);
    Etotal5=Etotal5+yf5(:,17+k);
    Etotal6=Etotal6+yf6(:,17+k);
    Etotal7=Etotal7+yf7(:,17+k);
    Etotal8=Etotal8+yf8(:,17+k);
    Etotal9=Etotal9+yf9(:,17+k);
    Istotal=Istotal+yf(:,34+k);
    Istotal2=Istotal2+yf2(:,34+k);
    Istotal3=Istotal3+yf3(:,34+k);
    Istotal4=Istotal4+yf4(:,34+k);
    Istotal5=Istotal5+yf5(:,34+k);
    Istotal6=Istotal6+yf6(:,34+k);
    Istotal7=Istotal7+yf7(:,34+k);
    Istotal8=Istotal8+yf8(:,34+k);
    Istotal9=Istotal9+yf9(:,34+k);
    Iatotal=Iatotal+yf(:,51+k);
    Iatotal2=Iatotal2+yf2(:,51+k);
    Iatotal3=Iatotal3+yf3(:,51+k);
    Iatotal4=Iatotal4+yf4(:,51+k);
    Iatotal5=Iatotal5+yf5(:,51+k);
    Iatotal6=Iatotal6+yf6(:,51+k);
    Iatotal7=Iatotal7+yf7(:,51+k);
    Iatotal8=Iatotal8+yf8(:,51+k);
    Iatotal9=Iatotal9+yf9(:,51+k);
    Htotal=Htotal+yf(:,68+k);
    Htotal2=Htotal2+yf2(:,68+k);
    Htotal3=Htotal3+yf3(:,68+k);
    Htotal4=Htotal4+yf4(:,68+k);
    Htotal5=Htotal5+yf5(:,68+k);
    Htotal6=Htotal6+yf6(:,68+k);
    Htotal7=Htotal7+yf7(:,68+k);
    Htotal8=Htotal8+yf8(:,68+k);
    Htotal9=Htotal9+yf9(:,68+k);
    Rtotal=Rtotal+yf(:,85+k);
    Rtotal2=Rtotal2+yf2(:,85+k);
    Rtotal3=Rtotal3+yf3(:,85+k);
    Rtotal4=Rtotal4+yf4(:,85+k);
    Rtotal5=Rtotal5+yf5(:,85+k);
    Rtotal6=Rtotal6+yf6(:,85+k);
    Rtotal7=Rtotal7+yf7(:,85+k);
    Rtotal8=Rtotal8+yf8(:,85+k);
    Rtotal9=Rtotal9+yf9(:,85+k);
    Dtotal=Dtotal+yf(:,102+k);
    Dtotal2=Dtotal2+yf2(:,102+k);
    Dtotal3=Dtotal3+yf3(:,102+k);
    Dtotal4=Dtotal4+yf4(:,102+k);
    Dtotal5=Dtotal5+yf5(:,102+k);
    Dtotal6=Dtotal6+yf6(:,102+k);
    Dtotal7=Dtotal7+yf7(:,102+k);
    Dtotal8=Dtotal8+yf8(:,102+k);
    Dtotal9=Dtotal9+yf9(:,102+k);
end

figure(1)

 

subplot(4,1,1)
plot(tf,Iatotal,'k','LineWidth',3)
hold on
plot(tf2,Iatotal2,'k','LineWidth',3)
plot(tf3,Iatotal3,'k','LineWidth',3)
plot(tf4,Iatotal4,'k','LineWidth',3)
plot(tf5,Iatotal5,'k','LineWidth',3)
plot(tf6,Iatotal6,'k','LineWidth',3)
plot(tf7,Iatotal7,'k','LineWidth',3)
plot(tf8,Iatotal8,'k','LineWidth',3)
plot(tf9,Iatotal9,'k','LineWidth',3)
ylabel('I_a','FontSize',16)
xticks([0 50 100 150 200 250 300 350 400 450])
xtickangle(45)

subplot(4,1,2)
plot(tf,Istotal,'k','LineWidth',3)
hold on
plot(tf2,Istotal2,'k','LineWidth',3)
plot(tf3,Istotal3,'k','LineWidth',3)
plot(tf4,Istotal4,'k','LineWidth',3)
plot(tf5,Istotal5,'k','LineWidth',3)
plot(tf6,Istotal6,'k','LineWidth',3)
plot(tf7,Istotal7,'k','LineWidth',3)
plot(tf8,Istotal8,'k','LineWidth',3)
plot(tf9,Istotal9,'k','LineWidth',3)
ylabel('I_s','FontSize',16)
xticks([0 50 100 150 200 250 300 350 400 450])
xtickangle(45)


subplot(4,1,3)
plot(tf,Htotal,'k','LineWidth',3)
hold on
plot(tf2,Htotal2,'k','LineWidth',3)
plot(tf3,Htotal3,'k','LineWidth',3)
plot(tf4,Htotal4,'k','LineWidth',3)
plot(tf5,Htotal5,'k','LineWidth',3)
plot(tf6,Htotal6,'k','LineWidth',3)
plot(tf7,Htotal7,'k','LineWidth',3)
plot(tf8,Htotal8,'k','LineWidth',3)
plot(tf9,Htotal9,'k','LineWidth',3)
ylabel('H','FontSize',16)
xticks([0 50 100 150 200 250 300 350 400 450])
xtickangle(45)

subplot(4,1,4)
plot(tf,Dtotal,'k','LineWidth',3)
hold on
plot(tf2,Dtotal2,'k','LineWidth',3)
plot(tf3,Dtotal3,'k','LineWidth',3)
plot(tf4,Dtotal4,'k','LineWidth',3)
plot(tf5,Dtotal5,'k','LineWidth',3)
plot(tf6,Dtotal6,'k','LineWidth',3)
plot(tf7,Dtotal7,'k','LineWidth',3)
plot(tf8,Dtotal8,'k','LineWidth',3)
plot(tf9,Dtotal9,'k','LineWidth',3)
ylabel('D','FontSize',16)
xticks([0 50 100 150 200 250 300 350 400 450])
xticklabels({'03/04','04/23','06/12','08/01','09/20','11/09','12/29','02/17','03/08','04/20'})
xtickangle(45)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p6,q6]=size(yf6);
IC7(1,:)=yf6(p6,:);

par.rho=repmat(0.0295,1,68);


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
for l=1:68
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M7(1,l)=max(eig(R2));
end
M7;


ft=linspace(83,150,68);
[tf7,yf7]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1):1:150],IC7);


[p7,q7]=size(yf7);
IC8(1,:)=yf7(p7,:);


par.rho=repmat(0.0295,1,310);


vw=1.0*ones(1,17); vo=1.0*ones(1,17); vh=1.0*ones(1,17);vs=1.0*ones(1,17);
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
M8=zeros(size(par.rho));
for l=1:310
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M8(1,l)=max(eig(R2));
end
M8;


ft=linspace(150,459,310);
[tf8,yf8]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1):1:459],IC8);

[p8,q8]=size(yf8);
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

YF=[yf; yf2; yf3; yf4; yf5; yf6; yf7; yf8];
Ro=[M M2 M3 M4 M5 M6 M7 M8];        
dlmwrite("lockdownJun1.csv",YF);                 
dlmwrite("lockdownRoJun1.csv",Ro);

Stotal=zeros(p,1);
Stotal2=zeros(p2,1);
Stotal3=zeros(p3,1);
Stotal4=zeros(p4,1);
Stotal5=zeros(p5,1);
Stotal6=zeros(p6,1);
Stotal7=zeros(p7,1);
Stotal8=zeros(p8,1);
Etotal=zeros(p,1);
Etotal2=zeros(p2,1);
Etotal3=zeros(p3,1);
Etotal4=zeros(p4,1);
Etotal5=zeros(p5,1);
Etotal6=zeros(p6,1);
Etotal7=zeros(p7,1);
Etotal8=zeros(p8,1);
Istotal=zeros(p,1);
Istotal2=zeros(p2,1);
Istotal3=zeros(p3,1);
Istotal4=zeros(p4,1);
Istotal5=zeros(p5,1);
Istotal5=zeros(p5,1);
Istotal6=zeros(p6,1);
Istotal7=zeros(p7,1);
Istotal8=zeros(p8,1);
Iatotal=zeros(p,1);
Iatotal2=zeros(p2,1);
Iatotal3=zeros(p3,1);
Iatotal4=zeros(p4,1);
Iatotal5=zeros(p5,1);
Iatotal6=zeros(p6,1);
Iatotal7=zeros(p7,1);
Iatotal8=zeros(p8,1);
Htotal=zeros(p,1);
Htotal2=zeros(p2,1);
Htotal3=zeros(p3,1);
Htotal4=zeros(p4,1);
Htotal5=zeros(p5,1);
Htotal5=zeros(p5,1);
Htotal6=zeros(p6,1);
Htotal7=zeros(p7,1);
Htotal8=zeros(p8,1);
Rtotal=zeros(p,1);
Rtotal2=zeros(p2,1);
Rtotal3=zeros(p3,1);
Rtotal4=zeros(p4,1);
Rtotal5=zeros(p5,1);
Rtotal6=zeros(p6,1);
Rtotal7=zeros(p7,1);
Rtotal8=zeros(p8,1);
Dtotal=zeros(p,1);
Dtotal2=zeros(p2,1);
Dtotal3=zeros(p3,1);
Dtotal4=zeros(p4,1);
Dtotal5=zeros(p5,1);
Dtotal5=zeros(p5,1);
Dtotal6=zeros(p6,1);
Dtotal7=zeros(p7,1);
Dtotal8=zeros(p8,1);
for k=1:17
    Stotal=Stotal+yf(:,0+k);
    Stotal2=Stotal2+yf2(:,0+k);
    Stotal3=Stotal3+yf3(:,0+k);
    Stotal4=Stotal4+yf4(:,0+k);
    Stotal5=Stotal5+yf5(:,0+k);
    Stotal6=Stotal6+yf6(:,0+k);
    Stotal7=Stotal7+yf7(:,0+k);
    Stotal8=Stotal8+yf8(:,0+k);
    Etotal=Etotal+yf(:,17+k);
    Etotal2=Etotal2+yf2(:,17+k);
    Etotal3=Etotal3+yf3(:,17+k);
    Etotal4=Etotal4+yf4(:,17+k);
    Etotal5=Etotal5+yf5(:,17+k);
    Etotal6=Etotal6+yf6(:,17+k);
    Etotal7=Etotal7+yf7(:,17+k);
    Etotal8=Etotal8+yf8(:,17+k);
    Istotal=Istotal+yf(:,34+k);
    Istotal2=Istotal2+yf2(:,34+k);
    Istotal3=Istotal3+yf3(:,34+k);
    Istotal4=Istotal4+yf4(:,34+k);
    Istotal5=Istotal5+yf5(:,34+k);
    Istotal6=Istotal6+yf6(:,34+k);
    Istotal7=Istotal7+yf7(:,34+k);
    Istotal8=Istotal8+yf8(:,34+k);
    Iatotal=Iatotal+yf(:,51+k);
    Iatotal2=Iatotal2+yf2(:,51+k);
    Iatotal3=Iatotal3+yf3(:,51+k);
    Iatotal4=Iatotal4+yf4(:,51+k);
    Iatotal5=Iatotal5+yf5(:,51+k);
    Iatotal6=Iatotal6+yf6(:,51+k);
    Iatotal7=Iatotal7+yf7(:,51+k);
    Iatotal8=Iatotal8+yf8(:,51+k);
    Htotal=Htotal+yf(:,68+k);
    Htotal2=Htotal2+yf2(:,68+k);
    Htotal3=Htotal3+yf3(:,68+k);
    Htotal4=Htotal4+yf4(:,68+k);
    Htotal5=Htotal5+yf5(:,68+k);
    Htotal6=Htotal6+yf6(:,68+k);
    Htotal7=Htotal7+yf7(:,68+k);
    Htotal8=Htotal8+yf8(:,68+k);
    Rtotal=Rtotal+yf(:,85+k);
    Rtotal2=Rtotal2+yf2(:,85+k);
    Rtotal3=Rtotal3+yf3(:,85+k);
    Rtotal4=Rtotal4+yf4(:,85+k);
    Rtotal5=Rtotal5+yf5(:,85+k);
    Rtotal6=Rtotal6+yf6(:,85+k);
    Rtotal7=Rtotal7+yf7(:,85+k);
    Rtotal8=Rtotal8+yf8(:,85+k);
    Dtotal=Dtotal+yf(:,102+k);
    Dtotal2=Dtotal2+yf2(:,102+k);
    Dtotal3=Dtotal3+yf3(:,102+k);
    Dtotal4=Dtotal4+yf4(:,102+k);
    Dtotal5=Dtotal5+yf5(:,102+k);
    Dtotal6=Dtotal6+yf6(:,102+k);
    Dtotal7=Dtotal7+yf7(:,102+k);
    Dtotal8=Dtotal8+yf8(:,102+k);
end

figure(1)


subplot(4,1,1)
hold on
plot(tf,Iatotal,'k','LineWidth',3)
hold on
plot(tf2,Iatotal2,'k','LineWidth',3)
plot(tf3,Iatotal3,'k','LineWidth',3)
plot(tf4,Iatotal4,'k','LineWidth',3)
plot(tf5,Iatotal5,'k','LineWidth',3)
plot(tf6,Iatotal6,'k','LineWidth',3)
plot(tf7,Iatotal7,'b','LineWidth',3)
plot(tf8,Iatotal8,'b','LineWidth',3)
ylabel('I_a','FontSize',16)
xticks([0 50 100 150 200 250 300 350 400 450])
set(gca,'Xticklabel',[])
xtickangle(45)
axis([0 450 0 300000])


subplot(4,1,2)
hold on
plot(tf,Istotal,'k','LineWidth',3)
hold on
plot(tf2,Istotal2,'k','LineWidth',3)
plot(tf3,Istotal3,'k','LineWidth',3)
plot(tf4,Istotal4,'k','LineWidth',3)
plot(tf5,Istotal5,'k','LineWidth',3)
plot(tf6,Istotal6,'k','LineWidth',3)
plot(tf7,Istotal7,'b','LineWidth',3)
plot(tf8,Istotal8,'b','LineWidth',3)
ylabel('I_s','FontSize',16)
xticks([0 50 100 150 200 250 300 350 400 450])
set(gca,'Xticklabel',[])
xtickangle(45)
axis([0 450 0 100000])


subplot(4,1,3)
hold on
plot(tf,Htotal,'k','LineWidth',3)
hold on
plot(tf2,Htotal2,'k','LineWidth',3)
plot(tf3,Htotal3,'k','LineWidth',3)
plot(tf4,Htotal4,'k','LineWidth',3)
plot(tf5,Htotal5,'k','LineWidth',3)
plot(tf6,Htotal6,'k','LineWidth',3)
plot(tf7,Htotal7,'b','LineWidth',3)
plot(tf8,Htotal8,'b','LineWidth',3)
ylabel('H','FontSize',16)
xticks([0 50 100 150 200 250 300 350 400 450])
set(gca,'Xticklabel',[])
xtickangle(45)
axis([0 450 0 10000])



subplot(4,1,4)
hold on
plot(tf,Dtotal,'k','LineWidth',3)
hold on
plot(tf2,Dtotal2,'k','LineWidth',3)
plot(tf3,Dtotal3,'k','LineWidth',3)
plot(tf4,Dtotal4,'k','LineWidth',3)
plot(tf5,Dtotal5,'k','LineWidth',3)
plot(tf6,Dtotal6,'k','LineWidth',3)
plot(tf7,Dtotal7,'b','LineWidth',3)
plot(tf8,Dtotal8,'b','LineWidth',3)
ylabel('D','FontSize',16)
xticks([0 50 100 150 200 250 300 350 400 450])
xticklabels({'2020-03-04','2020-04-23','2020-06-12','2020-08-01','2020-09-20','2020-11-09','2020-12-29','2021-02-17','2021-03-08','2021-04-20'})
xtickangle(45)
axis([0 450 0 40000])
