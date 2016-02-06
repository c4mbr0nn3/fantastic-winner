% Preconditioned Conjugate Gradient (PCG) (h=x0)
function [u,taureal,iter,Rres]=PCG(h,b,SYSMAT,IA,JA,N,NTERM,PREC,choice3)
if choice3==1
    DIAG=zeros(1,N);
    for i=1:N
        Spos=IA(i);
        DIAG(i)=1/SYSMAT(Spos);
        h(i)=DIAG(i)*b(i);
    end
end
% if choice2==1
[Ah]=matvec(N,IA,SYSMAT,h,JA);
r=b-Ah;
% end
u=h;
clear h
if choice3==1
    h=r.*DIAG;
elseif choice3==2
    h=avind(N,PREC,IA,JA,NTERM,r);
end
tau=norm(r)/norm(b);
tol=10e-15;
iter=0;
Rres=zeros(1,1000);
while (tau>tol) & (iter<1000)
    [Ah]=matvec(N,IA,SYSMAT,h,JA);
    alpha=(r*h.')/(h*Ah.');
    u=u+(alpha*h);
    r=r-(alpha*Ah);
    if choice3==1
        v=DIAG.*r;
    elseif choice3==2
        v=avind(N,PREC,IA,JA,NTERM,r);
    end
    beta=-((v*Ah.')/(h*Ah.'));
    h=v+(beta*h);
    iter=iter+1;
    tau=norm(r)/norm(b);
    Rres(1,iter)=tau;
end
Rres=Rres(Rres~=0);
u=u.';
% Calculate real residual
[Au]=matvec(N,IA,SYSMAT,u,JA);
taureal=norm(b-Au)/norm(b);
