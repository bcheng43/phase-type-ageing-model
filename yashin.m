rng(456)
lambda0=0.6;lambda=0.07;mu=0.4*10^-4;mu0=0.001;
yi=special(lambda0,lambda,mu0,mu,5000);
gridsize=500;
grid=max(yi)/gridsize;
y=[];
for ii=1:gridsize
    if sum(yi<=(grid*ii))-sum(yi<=(grid*(ii-1)))>0
   y=[y; [grid*ii, sum(yi<=(grid*ii))-sum(yi<=(grid*(ii-1)))]];
    end
end
profile on
tic
dx=[];
for ii=4:16
    dx=[dx;[estimatefun_fix(ii*25,y,0.999),ii*25]];
end
toc
profile viewer
[~,i]=min(dx(:,1));
x=dx(i,:);x(1)=[];m=x(end);%k=x(end-1);
n=m;
nn=1:n;
q1=x(1);qn=x(2);lam=x(3);tun=x(4);%qk=x(5);tun2=x(6);
if abs(tun)<10^-3
    qq=q1.^((n-nn)/(n-1)).*qn.^((nn-1)/(n-1));
else
    qq=(q1^tun*(n-nn)/(n-1)+qn^tun*(nn-1)/(n-1)).^(1/tun);
end
lambdaup=lam*ones(1,(n-1));
T=diag(-qq-[lambdaup,0]);T=T+diag(lambdaup,1);
sv=[];
sur_ph=[];
for i=1:length(y(:,1))
s=myexpmpdf(m,qq,lambdaup,y(i,1),10^-9);s=s(1);
sv=[sv,s];
sur=myexpm(m,qq,lambdaup,y(i,1),10^-9);sur=sur(1);
sur_ph=[sur_ph,sur];
end
%x=[lambda,mu]
specialsf=@(x) exp(-(lambda0+mu0)*x).*((lambda+mu)./(mu+lambda*exp(-(lambda+mu)*x))).^(lambda0/lambda);
s_data=specialsf(y(:,1));
specialhf=@(x) (mu0*mu+mu*lambda0+exp(-(lambda+mu)*x)*(mu0*lambda-mu*lambda0))./(mu+lambda*exp(-(lambda+mu)*x));  
specialf=@(x) exp(-(lambda0+mu0)*x).*((lambda+mu)./(mu+lambda*exp(-(lambda+mu)*x))).^(lambda0/lambda).*(mu0*mu+mu*lambda0+exp(-(lambda+mu)*x)*(mu0*lambda-mu*lambda0))./(mu+lambda*exp(-(lambda+mu)*x));
h_data=specialhf(y(:,1));
%emperical data
y_raw=[];
for i=1:size(y,1)
y_raw=[y_raw,y(i,1)*ones(1,y(i,2))];
end
[F,x]=ecdf(y_raw);
em_f=[0;diff(F)]./[1;diff(x)];

subplot(2,2,1)
plot(y(:,1),log(sur_ph),'black',y(:,1),log(s_data),'r','LineWidth',1.5)
hold on
% plot(x,log(1-F),'blue')
legend({'PHAM','Le Bras','Emperical'},'Location','northeast')
title('log(S(t))')

subplot(2,2,2)
plot(y(:,1),sv,'black',y(:,1),h_data.*s_data,'r','LineWidth',1.5)
hold on
% plot(x,em_f,'blue')
legend({'PHAM','Le Bras','Emperical'},'Location','northeast')
title('f(t)')

subplot(2,2,3)
plot(y(:,1),log(sv),'black',y(:,1),log(h_data.*s_data),'r','LineWidth',1.5)
hold on
% plot(x,log(em_f),'blue')
legend({'PHAM','Le Bras','Emperical'},'Location','northeast')
title('log(f(t))')

subplot(2,2,4)
plot(y(:,1),sv./sur_ph,'black',y(:,1),h_data,'r','LineWidth',1.5)
hold on
% plot(x,em_f./(1-F),'blue')
legend({'PHAM','Le Bras','Emperical'},'Location','northeast')
title('hazard')


