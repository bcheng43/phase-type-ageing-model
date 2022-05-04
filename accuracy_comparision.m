format long
%linear hi 
b=0.015;
c=0.01;
m=50;
nn=1:m;
qq_real=b+c*nn;
f_theoretical=[];
f_built=[];
f_my=[];
for t=1:40
lambda=1.6; 
%theoretical
p=zeros(m,1);
p(1)=exp(-(qq_real(1)+lambda)*t);
con=lambda/c*(1-exp(-c*t));
for k=2:(m-1)
p(k)=p(k-1)*con/(k-1);
end
fun=@(s) exp(-(lambda+qq_real(1)-qq_real(m))*s).*(lambda/c*(1-exp(-c*s))).^(m-2);
int=integral(fun,0,t,'RelTol',0,'AbsTol',1e-12);
p(m)=int/factorial(m-2)*exp(-qq_real(m)*t)*lambda;
f_theo=qq_real*p;
f_theoretical=[f_theoretical,f_theo];
%bulit in
T=diag(-lambda-qq_real);
T=T+diag(lambda*ones(1,m-1),1);
f=expm(T*t)*qq_real'; f=f(1);
f_built=[f_built,f];
%
err_d=10^-10;
y=t;
%%%%%%%%%%%%%%%%
lambda=lambda*ones(1,m-1);
maxx=1.001*max(qq_real+[lambda,0]);
%find truncation N 
N=0;
sum_N=0;
P=zeros(1,m);
P(1)=1-(lambda(1)+qq_real(1))/maxx;
P(2)=lambda(1)/maxx;
kk=poisspdf(0,maxx*y(:,1));
ex=zeros(size(y,1),m);
ex(:,1)=kk;
 while (sum_N<=(1-err_d))
kk=poisspdf(N+1,maxx*y(:,1));
sum_N=sum_N+poisspdf(N,maxx*y(end,1));
N=N+1;
ex=ex+kk*P;
%modify recursive formula
PP=[0,P];
PP(end)=[];
PP=[0,lambda]/maxx.*PP+(1-([lambda,0]+qq_real)/maxx).*P;
% PP=[0,lambda]/maxx.*PP+(1-([lambda,x(3)]+qq)/maxx).*P;
P=PP;
 end
like=ex*qq_real';
f_my=[f_my,like];
end

h = figure;
plot(f_theoretical-f_built)
hold on
plot(f_theoretical-f_my)
xlabel('t')
ylabel('difference between two pdf')
legend({'f1-f2','f1-f3'},'Location','northwest')
title('f1 theoretical pdf; f2 pdf using expm; f3 pdf caculated from our algorithm')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'filename','-dpdf','-r0')
