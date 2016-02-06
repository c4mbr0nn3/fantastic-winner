function [h]=CR(b,w,N,IA,SYSMAT,maxiter,JA,PREC,NTERM)
iter=0;
h=w;
[Ah]=matvec(N,IA,SYSMAT,h,JA);
v=b-Ah;
while iter<maxiter
    [w]=avind(N,PREC,IA,JA,NTERM,v);
    h=h+w;
    [Ah]=matvec(N,IA,SYSMAT,h,JA);
    v=b-Ah;
    iter=iter+1;
end