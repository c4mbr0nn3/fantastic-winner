function [h]=CR(b,N,IA,SYSMAT,maxiter,JA,PREC,NTERM)
h=avind(N,PREC,IA,JA,NTERM,b); %(h=x0)
for i=1:maxiter
    r=b-matvec(N,IA,SYSMAT,h,JA);
    h=h+avind(N,PREC,IA,JA,NTERM,r);
    
end