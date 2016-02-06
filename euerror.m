% Error function
function [error,ansol,perror]=euerror(N,NCOORD,areanod,u)
afun=@(x,y)(x.^2+y.^2-(x.^2.*y.^2)-1);
serror=0;
ansol=zeros(1,N);
perror=zeros(1,N);
for i=1:N
    Xnode=NCOORD(i,1);
    Ynode=NCOORD(i,2);
    serror=serror+(afun(Xnode,Ynode)-u(i))^2*areanod(i);
    ansol(i)=afun(Xnode,Ynode);
    perror(i)=abs(u(i)-afun(Xnode,Ynode));
end
error=sqrt(serror);