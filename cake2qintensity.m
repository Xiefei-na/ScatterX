

function [q,inten] = cake2qintensity(A,cenx,ceny,phi1,phi2,rmax,lambda,sampletodetector,pixelsize)

A(A<1)=NaN;
[lx,la]=size(A);
interpolatedA=@(xx,yy) interp2(1:la,1:lx,A,xx,yy);
r = 1:rmax;
theta=phi1:phi2;
theta=theta+180;theta_rad = theta*pi/180;
x = r'*cos(theta_rad)+cenx;
y = r'*sin(theta_rad)+ceny;
Si=interpolatedA(x,y);Si=Si';

inten = nanmean(Si, 1);
q = 4*pi*sin(atan((1:rmax)*pixelsize/sampletodetector)/2)/lambda;

end