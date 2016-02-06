% Calculate known term f without b.c.
function [b,SYSMAT]=bcond(N,NCOORD,areanod,NDnum,NNnum,NNbound,NDbound,IA,SYSMAT)
b0=zeros(1,N);
% syms f(x,y)
% prompt0='Equazione f: ';
% f(x,y)=input(prompt0);
fun=@(x,y)(- 2.*x.^2 - 2.*y.^2 + 4);
for k=1:N
    Xnode=NCOORD(k,1);
    Ynode=NCOORD(k,2);
%     b0(k)=f(Xnode,Ynode)*areanod(k);
    b0(k)=fun(Xnode,Ynode)*areanod(k);
    clear('Xnode','Ynode')
end

% Apply Neumann b.condition
neucorr=zeros(1,N);
l=4/((NDnum/2)-1); 
fun2=@(y)(y.^2-1);
for k=1:NNnum
node=NNbound(k); 
Ynode=NCOORD(node,2); 
neucorr(node)=fun2(Ynode)*l;
clear('node','Ynode')
end
b=b0+neucorr;

% Apply Dirichlet b.condition
for k=1:NDnum
node=NDbound(k);
b(node)=0;
ind=IA(node);
SYSMAT(ind)=10e15;
clear('node','ind');
end
