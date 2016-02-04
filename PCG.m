% Preconditioned Conjugate Gradient (PCG) (h=x0)
function [x,tau,taureal,iter]=PCG(h,b,SYSMAT,IA,JA,N,NTERM,PREC)
[Ah]=matvec(N,IA,SYSMAT,h,JA);
x=h;
r=b-Ah;
v=r;
[w]=avind(N,PREC,IA,JA,NTERM,v); % posso usare direttamente p??
h=w;
tau=norm(r)/norm(b);
tol=10e-15;
iter=0;
while tau>tol
    if iter>1000
        break
    else
        [Ah]=matvec(N,IA,SYSMAT,h,JA);
        alpha=(r*h.')/(h*Ah.');
        x=x+(alpha*h);
        r=r-(alpha*Ah);
        v=r;
        [w]=avind(N,PREC,IA,JA,NTERM,v);
        v=w;
        beta=-((v*Ah.')/(h*Ah.'));
        h=v+(beta*h);
        iter=iter+1;
        tau(iter)=norm(r)/norm(b);
    end
end

% Calculate real residual
h=x;
[Ah]=matvec(N,IA,SYSMAT,h,JA);
taureal=norm(b-Ah)/norm(b);
