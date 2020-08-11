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

par.rho=repmat(0.0295,1,11);


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
for l=1:11
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M7(1,l)=max(eig(R2));
end
M7;


ft=linspace(83,93,11);
[tf7,yf7]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1):1:93],IC7);


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

[p7,q7]=size(yf7);
IC8(1,:)=yf7(p7,:);


par.rho=repmat(0.0295,1,11);


vw=0.3*ones(1,17); vo=0.3*ones(1,17); vh=1.0*ones(1,17);vs=0.0*ones(1,17);
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
for l=1:11
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M8(1,l)=max(eig(R2));
end
M8;


ft=linspace(93,103,11);
[tf8,yf8]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1):1:103],IC8);

[p8,q8]=size(yf8);
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
IC9(1,:)=yf8(p8,:);

par.rho=repmat(0.0295,1,11);

vw=0.4*ones(1,17); vo=0.4*ones(1,17); vh=1.0*ones(1,17);
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
M9=zeros(size(par.rho));
for l=1:11
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M9(1,l)=max(eig(R2));
end
M9;

ft=linspace(103,113,11);
[tf9,yf9]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1)+(p8-1):1:113],IC9);


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



[p9,q9]=size(yf9);
IC10(1,:)=yf9(p9,:);

par.rho=repmat(0.0295,1,11);

vw=0.5*ones(1,17); vo=0.5*ones(1,17); vh=1.0*ones(1,17);
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
M10=zeros(size(par.rho));
for l=1:11
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M10(1,l)=max(eig(R2));
end
M10;

ft=linspace(113,123,11);
[tf10,yf10]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1)+(p8-1)+(p9-1):1:123],IC10);




%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



[p10,q10]=size(yf10);
IC11(1,:)=yf10(p10,:);

par.rho=repmat(0.0295,1,11);

vw=0.6*ones(1,17); vo=0.6*ones(1,17); vh=1.0*ones(1,17);
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
M11=zeros(size(par.rho));
for l=1:11
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M11(1,l)=max(eig(R2));
end
M11;

ft=linspace(123,133,11);
[tf11,yf11]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1)+(p8-1)+(p9-1)+(p10-1):1:133],IC11);


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



[p11,q11]=size(yf11);
IC12(1,:)=yf11(p11,:);

par.rho=repmat(0.0295,1,11);

vw=0.7*ones(1,17); vo=0.7*ones(1,17); vh=1.0*ones(1,17);
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
M12=zeros(size(par.rho));
for l=1:11
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M12(1,l)=max(eig(R2));
end
M12;

ft=linspace(133,143,11);
[tf12,yf12]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1)+(p8-1)+(p9-1)+(p10-1)+(p11-1):1:143],IC12);


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p12,q12]=size(yf12);
IC13(1,:)=yf12(p12,:);

par.rho=repmat(0.0295,1,11);

vw=0.8*ones(1,17); vo=0.8*ones(1,17); vh=1.0*ones(1,17);
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
M13=zeros(size(par.rho));
for l=1:11
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M13(1,l)=max(eig(R2));
end
M13;

ft=linspace(143,153,11);
[tf13,yf13]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1)+(p8-1)+(p9-1)+(p10-1)+(p11-1)+(p12-1):1:153],IC13);



%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p13,q13]=size(yf13);
IC14(1,:)=yf13(p13,:);

par.rho=repmat(0.0295,1,11);

vw=0.9*ones(1,17); vo=0.9*ones(1,17); vh=1.0*ones(1,17);
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
M14=zeros(size(par.rho));
for l=1:11
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M14(1,l)=max(eig(R2));
end
M14;

ft=linspace(153,163,11);
[tf14,yf14]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1)+(p8-1)+(p9-1)+(p10-1)+(p11-1)+(p12-1)+(p13-1):1:163],IC14);


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p14,q14]=size(yf14);
IC15(1,:)=yf14(p14,:);

par.rho=repmat(0.0295,1,19);

vw=1.0*ones(1,17); vo=1.0*ones(1,17); vh=1.0*ones(1,17);
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
M15=zeros(size(par.rho));
for l=1:19
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M15(1,l)=max(eig(R2));
end
M15;


ft=linspace(163,181,19);
[tf15,yf15]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1)+(p8-1)+(p9-1)+(p10-1)+(p11-1)+(p12-1)+(p13-1)+(p14-1):1:181],IC15);


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


[p15,q15]=size(yf15);
IC16(1,:)=yf15(p15,:);

par.rho=repmat(0.0295,1,270);

vw=1.0*ones(1,17); vo=1.0*ones(1,17); vh=1.0*ones(1,17); vs=1.0*ones(1,17);
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
M16=zeros(size(par.rho));
for l=1:270
R2=par.rho(l)/((par.sigmaa+par.sigmas)/2)*Crec2;
M16(1,l)=max(eig(R2));
end
M16;

ft=linspace(181,450,270);
[tf16,yf16]=ode45(@(t,y)covidnj(t,y,par,age,Contact,ft),[(p-1)+(p2-1)+(p3-1)+(p4-1)+(p5-1)+(p6-1)+(p7-1)+(p8-1)+(p9-1)+(p10-1)+(p11-1)+(p12-1)+(p13-1)+(p14-1)+(p15-1):1:450],IC16);

[p16,q16]=size(yf16);


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

YF=[yf; yf2; yf3; yf4; yf5; yf6; yf7; yf8; yf9; yf10; yf11; yf12; yf13; yf14; yf15; yf16];
Ro=[M M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12 M13 M14 M15 M16]; 
dlmwrite("lockdowncontinuedrelaxation2.csv",YF);                 
dlmwrite("lockdownRocontinuedrelaxation2.csv",Ro);



Stotal=zeros(p,1);
Stotal2=zeros(p2,1);
Stotal3=zeros(p3,1);
Stotal4=zeros(p4,1);
Stotal5=zeros(p5,1);
Stotal6=zeros(p6,1);
Stotal7=zeros(p7,1);
Etotal=zeros(p,1);
Etotal2=zeros(p2,1);
Etotal3=zeros(p3,1);
Etotal4=zeros(p4,1);
Etotal5=zeros(p5,1);
Etotal6=zeros(p6,1);
Etotal7=zeros(p7,1);
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
Istotal10=zeros(p10,1);
Istotal11=zeros(p11,1);
Istotal12=zeros(p12,1);
Istotal13=zeros(p13,1);
Istotal14=zeros(p14,1);
Istotal15=zeros(p15,1);
Istotal16=zeros(p16,1);
Iatotal=zeros(p,1);
Iatotal2=zeros(p2,1);
Iatotal3=zeros(p3,1);
Iatotal4=zeros(p4,1);
Iatotal5=zeros(p5,1);
Iatotal6=zeros(p6,1);
Iatotal7=zeros(p7,1);
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
Htotal10=zeros(p10,1);
Htotal11=zeros(p11,1);
Htotal12=zeros(p12,1);
Htotal13=zeros(p13,1);
Htotal14=zeros(p14,1);
Htotal15=zeros(p15,1);
Htotal16=zeros(p16,1);
Rtotal=zeros(p,1);
Rtotal2=zeros(p2,1);
Rtotal3=zeros(p3,1);
Rtotal4=zeros(p4,1);
Rtotal5=zeros(p5,1);
Rtotal6=zeros(p6,1);
Rtotal7=zeros(p7,1);
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
Dtotal10=zeros(p10,1);
Dtotal11=zeros(p11,1);
Dtotal12=zeros(p12,1);
Dtotal13=zeros(p13,1);
Dtotal14=zeros(p14,1);
Dtotal15=zeros(p15,1);
Dtotal16=zeros(p16,1);
for k=1:17
    Stotal=Stotal+yf(:,0+k);
    Stotal2=Stotal2+yf2(:,0+k);
    Stotal3=Stotal3+yf3(:,0+k);
    Stotal4=Stotal4+yf4(:,0+k);
    Stotal5=Stotal5+yf5(:,0+k);
    Stotal6=Stotal6+yf6(:,0+k);
    Stotal7=Stotal7+yf7(:,0+k);
    Etotal=Etotal+yf(:,17+k);
    Etotal2=Etotal2+yf2(:,17+k);
    Etotal3=Etotal3+yf3(:,17+k);
    Etotal4=Etotal4+yf4(:,17+k);
    Etotal5=Stotal5+yf5(:,17+k);
    Etotal6=Stotal6+yf6(:,17+k);
    Etotal7=Stotal7+yf7(:,17+k);
    Istotal=Istotal+yf(:,34+k);
    Istotal2=Istotal2+yf2(:,34+k);
    Istotal3=Istotal3+yf3(:,34+k);
    Istotal4=Istotal4+yf4(:,34+k);
    Istotal5=Istotal5+yf5(:,34+k);
    Istotal6=Istotal6+yf6(:,34+k);
    Istotal7=Istotal7+yf7(:,34+k);
    Istotal8=Istotal8+yf8(:,34+k);
    Istotal9=Istotal9+yf9(:,34+k);
    Istotal10=Istotal10+yf10(:,34+k);
    Istotal11=Istotal11+yf11(:,34+k);
    Istotal12=Istotal12+yf12(:,34+k);
    Istotal13=Istotal13+yf13(:,34+k);
    Istotal14=Istotal14+yf14(:,34+k);
    Istotal15=Istotal15+yf15(:,34+k);
    Istotal16=Istotal16+yf16(:,34+k);
    Iatotal=Iatotal+yf(:,51+k);
    Iatotal2=Iatotal2+yf2(:,51+k);
    Iatotal3=Iatotal3+yf3(:,51+k);
    Iatotal4=Iatotal4+yf4(:,51+k);
    Iatotal5=Iatotal5+yf5(:,51+k);
    Iatotal6=Iatotal6+yf6(:,51+k);
    Iatotal7=Iatotal7+yf7(:,51+k);
    Htotal=Htotal+yf(:,68+k);
    Htotal2=Htotal2+yf2(:,68+k);
    Htotal3=Htotal3+yf3(:,68+k);
    Htotal4=Htotal4+yf4(:,68+k);
    Htotal5=Htotal5+yf5(:,68+k);
    Htotal6=Htotal6+yf6(:,68+k);
    Htotal7=Htotal7+yf7(:,68+k);
    Htotal8=Htotal8+yf8(:,68+k);
    Htotal9=Htotal9+yf9(:,68+k);
    Htotal10=Htotal10+yf10(:,68+k);
    Htotal11=Htotal11+yf11(:,68+k);
    Htotal12=Htotal12+yf12(:,68+k);
    Htotal13=Htotal13+yf13(:,68+k);
    Htotal14=Htotal14+yf14(:,68+k);
    Htotal15=Htotal15+yf15(:,68+k);
    Htotal16=Htotal16+yf16(:,68+k);
    Rtotal=Rtotal+yf(:,85+k);
    Rtotal2=Rtotal2+yf2(:,85+k);
    Rtotal3=Rtotal3+yf3(:,85+k);
    Rtotal4=Rtotal4+yf4(:,85+k);
    Rtotal5=Rtotal5+yf5(:,85+k);
    Rtotal6=Rtotal6+yf6(:,85+k);
    Rtotal7=Rtotal7+yf7(:,85+k);
    Dtotal=Dtotal+yf(:,102+k);
    Dtotal2=Dtotal2+yf2(:,102+k);
    Dtotal3=Dtotal3+yf3(:,102+k);
    Dtotal4=Dtotal4+yf4(:,102+k);
    Dtotal5=Dtotal5+yf5(:,102+k);
    Dtotal6=Dtotal6+yf6(:,102+k);
    Dtotal7=Dtotal7+yf7(:,102+k);
    Dtotal8=Dtotal8+yf8(:,102+k);
    Dtotal9=Dtotal9+yf9(:,102+k);
    Dtotal10=Dtotal10+yf10(:,102+k);
    Dtotal11=Dtotal11+yf11(:,102+k);
    Dtotal12=Dtotal12+yf12(:,102+k);
    Dtotal13=Dtotal13+yf13(:,102+k);
    Dtotal14=Dtotal14+yf14(:,102+k);
    Dtotal15=Dtotal15+yf15(:,102+k);
    Dtotal16=Dtotal16+yf16(:,102+k);
end



 figure(1)
 plot(tf,Istotal,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
hold on
plot(tf2,Istotal2,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf3,Istotal3,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf4,Istotal4,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf5,Istotal5,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf6,Istotal6,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf7,Istotal7,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf8,Istotal8,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf9,Istotal9,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf10,Istotal10,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf11,Istotal11,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf12,Istotal12,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf13,Istotal13,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf14,Istotal14,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
plot(tf15,Istotal15,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)
p1=plot(tf16,Istotal16,'Color',[0.3466 0.5360 0.6906],'LineWidth',3)

plot(tf,Htotal,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf2,Htotal2,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf3,Htotal3,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf4,Htotal4,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf5,Htotal5,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf6,Htotal6,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf7,Htotal7,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf8,Htotal8,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf9,Htotal9,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf10,Htotal10,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf11,Htotal11,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf12,Htotal12,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf13,Htotal13,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf14,Htotal14,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
plot(tf15,Htotal15,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)
p2=plot(tf16,Htotal16,'Color',[0.9152 0.2815 0.2878],'LineWidth',3)

plot(tf,Dtotal,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf2,Dtotal2,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf3,Dtotal3,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf4,Dtotal4,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf5,Dtotal5,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf6,Dtotal6,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf7,Dtotal7,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf8,Dtotal8,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf9,Dtotal9,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf10,Dtotal10,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf11,Dtotal11,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf12,Dtotal12,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf13,Dtotal13,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf14,Dtotal14,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
plot(tf15,Dtotal15,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
p3=plot(tf16,Dtotal16,'Color',[0.4415 0.7490 0.4321],'LineWidth',3)
legend([p1 p2 p3],'Symptomatic Infectious','Hospitalised','Cumulative Deaths','Location','NorthEast') 
xticks([0 50 100 150 200 250 300 350 400 450])
xticklabels({'2020-03-04','2020-04-23','2020-06-12','2020-08-01','2020-09-20','2020-11-09','2020-12-29','2021-02-17','2021-03-08','2021-04-20'})
xtickangle(45)
