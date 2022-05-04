% m total number of state; first column of y is data observartion, second
% column is number of observation; lquantile 1-alpha in TVaR
function [estimatefun]=estimatefun_fix(m,y,lquantile)
% naive try
% xr=zeros(50,5);
xr=zeros(10,8);
% xr=zeros(20,12);
y_raw=[];
for ii=1:size(y,1)
    y_raw=[y_raw,y(ii,1)*ones(1,y(ii,2))];
end
psi=mean(y_raw(y_raw>=quantile(y_raw,lquantile)));
% objfun2=@(x) objfun1(m,x(1),x(2),m/median(y_raw(y_raw>=quantile(y_raw,0.999))),x(3),x(4),x(5),k,y);
% input (total number of state, h1,hm,lambda,s,data,tolerance)
% objfun2=@(x) objfun1(m,x(1),x(2),m/mean(y_raw(y_raw>=quantile(y_raw,0.999))),x(3),y,10^-10);
% objfun2=@(x) objfun1(m,x(1),x(2),m/mean(y_raw(y_raw>=quantile(y_raw,lquantile))),x(3),y,10^-10);
objfun2=@(x) -(y(:,2)'*log(PTAM(y,m,x(1),x(2),m/psi,x(3),10^-10)));
% lb=[0,0,-100];ub=[40,40,100];A=[];B=[];
 lb=[0,0,-100];ub=[100,100,100];A=[];B=[];
x0_set=rand(size(xr,1),3);
x0_set(:,3)=x0_set(:,3)*2-1;
tic
parfor ii=1:size(x0_set,1)
x0=x0_set(ii,:);
% A=[1,m,0;-1,-m,0;1,1,0;-1,-1,0];B=[1;1;1;1]; ,'UseParallel',true
[x, fval,~,~,~,~,hessian]=fmincon(objfun2,x0,A,B,[],[],lb,ub,[],optimoptions('fmincon','MaxIter',10000,'MaxFunEvals',10000,'TolFun',1e-5,'TolX',1e-5));
% xr(ii,:)=[fval,x(1:2),m/mean(y_raw(y_raw>=quantile(y_raw,0.999))),x(3),x(4),x(5),eig(hessian)'];
% xr(ii,:)=[fval,x(1:2),m/mean(y_raw(y_raw>=quantile(y_raw,0.999))),x(3),eig(hessian)'];
xr(ii,:)=[fval,x(1:2),m/mean(y_raw(y_raw>=quantile(y_raw,lquantile))),x(3),eig(hessian)'];
% xr(ii,:)=[fval,x(1:2),m/mean(y_raw(y_raw>=quantile(y_raw,0.99))),x(3),eig(hessian)'];
end
toc
[~,i]=min(xr(:,1));

estimatefun=xr(i,:);

end