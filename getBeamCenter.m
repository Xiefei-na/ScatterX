function [cenx,ceny] = getBeamCenter(p)
[m,n]=find(   p  ==  round(max(p(:))/50)   );
X=[m,n];
n=length(X(:,1));
y=ones(n,1);
b=[rand(1)*1000 rand(1) rand(1)];
fun=inline('X(:,1).^2+X(:,2).^2+a(1)*X(:,1)+a(2)*X(:,2)+a(3)','a','X');
[a,~,~] = nlinfit(X,y,fun,b);
ceny =-a(1)/2;cenx =-a(2)/2;
end
