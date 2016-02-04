function [Ah]=matvec(N,IA,SYSMAT,h,JA)
Ah=zeros(1,N);
for i=1:N
    k=IA(i);
    Ah(i)=Ah(i)+SYSMAT(k)*h(i);
    for k=IA(i)+1:IA(i+1)-1
        j=JA(k);
        Ah(i)=Ah(i)+SYSMAT(k)*h(j);
        Ah(j)=Ah(j)+SYSMAT(k)*h(i);
    end
end