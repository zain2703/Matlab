function [lambda,v,Iter]=eig_p(A,tol)
    x0=ones(size(A,1),1);       %Inital guess for eigen value
    J=size(A,1)
    y1=tol+1;                   %variable to initlize values
    diff=100;                   %variable to initlize values 
    Iter=0;                     %Iteration values
while diff>tol
    y=x0/norm(x0);         %Normalizing first eigen vector
    x0=A*y;                %Next eigen vector
    y0=y'*x0;              %EigenValue
    diff=abs(y1-y0);       
    y1=y0;                 %Swapping eigen value                 
    Iter=Iter+1;
end
    lambda=y0;
    v=y;
end