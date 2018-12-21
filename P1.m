lambda=zeros(5,1);
Iter1=zeros(5,1);
Iter2=zeros(5,1);
n=linspace(20,100,5);
for i=1:5
    h=1/n(i);               %   Ratio at n Linespace
    K1=2*ones(1,n(i));      %   diagnol with 2
    K1(1,n(i))=1;           %   diagnol 2 with 1 at end
    K2=-1*ones(1,n(i)-1);   %   -1 diagnol
    M1=2*K1;                %   2*k1 to make 4 diagnol for m matrix 
    M1(1,n(i))=2;           %   diagnol with 2 at the end
    M2=-1*K2;               %   diagnol with 1 
    K=(diag(K1,0)+diag(K2,1)+diag(K2,-1))/h;     % K matrix
    M=h*(diag(M1,0)+diag(M2,1)+diag(M2,-1))/6;   % M matrix
    A=M\K;                                       % Matrix A
    [lambdaSmall(i),Vsmall,Iter1(i)]=eiginv_power(A,10e-4);
    [lambdaLarge(i),Vlarge,Iter2(i)]=eig_p(A,10e-4);
    MatEig(i)=max(eig(A));
end
   figure (1);clf
   plot(Vlarge,'r-')
   ylabel('Eigenvector n=100')
   xlabel('Element')
   hold on 
   plot(Vsmall,'b-')
   
   figure (2);clf
   plot(n,lambdaLarge,'r-o')
   ylabel('\Lambda')
   xlabel('n')
   hold on
   plot(n,MatEig,'b-x')