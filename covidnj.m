function dydt = covidnj(t,y,par,age,Contact,ft)


dydt=zeros(length(y),1);

par.rho = interp1(ft,par.rho,t);

lambda(age.S)=(par.rho*sum(Contact.Cf(age.S,1:17)./age.N(1,1:17).*y(age.Is)) + par.rho*sum(Contact.Cf(age.S,1:17)./age.N(1,1:17).*y(age.Ia)));

dydt(age.S)=-lambda(age.S).*y(age.S)';
dydt(age.E)=lambda(age.S).*y(age.S)'-par.gamma*y(age.E)';
dydt(age.Is)=par.gamma*age.f.*y(age.E)'-par.sigmas*y(age.Is)';
dydt(age.Ia)=par.gamma*(1-age.f).*y(age.E)'-par.sigmaa*y(age.Ia)';
dydt(age.H)=par.sigmas*age.h.*y(age.Is)'-par.alphah*y(age.H)';
dydt(age.R)=par.sigmas*(1-age.h-age.d).*y(age.Is)'+par.sigmaa*y(age.Ia)'+par.alphah*(1-age.dh).*y(age.H)';
dydt(age.D)=par.sigmas*age.d.*y(age.Is)'+par.alphah*age.dh.*y(age.H)';

