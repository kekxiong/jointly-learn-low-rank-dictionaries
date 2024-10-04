function U=UIQI(x,y)

[rows,columns,bands]=size(x);

x_=mean2(x);y_=mean2(y);
Sigma_x2=sum(sum(sum((x-x_))))/(rows*columns*bands-1);
Sigma_y2=sum(sum(sum((y-y_))))/(rows*columns*bands-1);
Sigma_xy=sum(sum(sum((x-x_).*(y-y_))))/(rows*columns*bands-1);

U1=Sigma_xy/(sqrt(Sigma_x2)*sqrt(Sigma_y2));
U2=(2*x_*y_)/(x_^2+y_^2);
U3=(2*sqrt(Sigma_x2)*sqrt(Sigma_y2))/(Sigma_x2+Sigma_y2);
U=U1*U2*U3;