% generate from survaival funtion. A duality in aging paper

function [special]=special(lambda0,lambda,mu0,mu,num)
% u=rand(1,num);
% rate=lambda+mu;
% x=log(rate-lambda*u)-log(mu*u);x=x/rate;
% special=x;

specialpdf=@(x) exp(-(lambda0+mu0)*x).*((lambda+mu)./(mu+lambda*exp(-(lambda+mu)*x))).^(lambda0/lambda);
age=1:20000;age=age*0.01;
survival=specialpdf(age);
u=rand(1,num);u=sort(u,'descend');
special=[];
for ii=1:num
    location=sum(survival>u(ii));
    special=[special;location*0.01];
end

end