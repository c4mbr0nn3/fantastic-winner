function [w]=avind(N,PREC,IA,JA,NTERM,v)
s=zeros(1,N);
z=zeros(1,N);
z(1)=v(1)/PREC(1);
for i=2:N
    k=IA(i);
    for m=IA(i-1)+1:IA(i)-1
        j=JA(m);
        s(j)=s(j)+PREC(m)*z(i-1);
    end
    z(i)=(v(i)-s(i))/PREC(k);
end
w=zeros(1,N);
% clear s
% s=zeros(1,N-1);
w(N)=z(N)/PREC(NTERM);
for i=N-1:-1:1
    s(i)=0;
    k=IA(i);
    for m=IA(i)+1:IA(i+1)-1
        j=JA(m);
        s(i)=s(i)+PREC(m)*w(j);
    end
    w(i)=(z(i)-s(i))/PREC(k);
end