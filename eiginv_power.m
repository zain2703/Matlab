function [lambda,v,iter] = eiginv_power(A,tol)
[rows,~]=size(A);
[L, U, P]= lu(A);        %LU fACT
x1=ones(rows,1);         %Guessing first eigen value
diff=tol+1;           
iter=0;
while diff>tol
    d=L\(P*x1);          
    x1=U\d;              
    %Next eigen vector
    x2=x1/norm(x1);      %Normalizing second eigen vector   
    A1 = A*x2;              %Next eigen vector
    m=(x2'*A1/(x2'*x2));    %Rayleigh quotient
    diff=norm(A1-m*x2);     
    iter=iter+1;
lambda = m;
disp(m)
end
%I=A1;           Mistakes
v=x2;
 disp('Smallest eigenvalue=');
 disp(lambda);
