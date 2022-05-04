function [obj]=PTAM(y,n,q1,qn,lam,tun,err_d)
%number of state m=n
nn=(2:n-1)'; 
if abs(tun)<10^-3
    qq=q1.^((n-nn)/(n-1)).*qn.^((nn-1)/(n-1));
else
    qq=(q1^tun*(n-nn)/(n-1)+qn^tun*(nn-1)/(n-1)).^(1/tun);
end
qq=[q1;qq;qn];
lambda=lam*ones((n-1),1);
%%%%%%%%%%%%%%%%
maxx=1.01*max(qq+[lambda;0]);
P=zeros(n,1);
P(1)=1-(lambda(1)+qq(1))/maxx;
P(2)=lambda(1)/maxx;
ex=zeros(size(y,1),n);
ex(:,1)=poisspdf(0,maxx*y(:,1));
%%%%%% uniformization and vectorization
%%truncation relative difference is less than tolerance
%weight for calculation PP
w1=[0;lambda]/maxx; w2=1-([lambda;0]+qq)/maxx;
%find truncation N
N=poissinv(1-err_d,maxx*y(end,1));
% N=sum(poisscdf((1:10*max(floor(maxx*y(end,1)),1))',maxx*y(end,1))<1-err_d);
% %in case N is not large enough
% if N==10*max(floor(maxx*y(end,1)),1)
% N=sum((1-poisscdf(1:100*max(floor(maxx*y(end,1)),1),maxx*y(end,1)))>err_d);    
% end
 for N_i=1:N
ex=ex+poisspdf(N_i,maxx*y(:,1))*P';
P=w1.*[0;P(1:end-1)]+w2.*P;
 end
obj=ex*qq;
end