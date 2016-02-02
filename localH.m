% Local Stiffness Matrix function
function [localH]=localH(TRIANG,NCOORD,NT)
localH=zeros(3,3,NT);
for e=1:NT
    node=zeros(1,3);
    Xnode=zeros(1,3);
    Ynode=zeros(1,3);
    for k=1:3
    node(k)=TRIANG(k,e);
    Xnode(k)=NCOORD(node(k),1).';
    Ynode(k)=NCOORD(node(k),2).';
    end
    b(1)=Ynode(2)-Ynode(3);
    b(2)=Ynode(3)-Ynode(1);
    b(3)=Ynode(2)-Ynode(1);
    c(1)=Xnode(3)-Xnode(2);
    c(2)=Xnode(1)-Xnode(3);
    c(3)=Xnode(2)-Xnode(1);
    D=det([ones(1,3);Xnode;Ynode].')/2;
    matB=zeros(3,3);
    matC=zeros(3,3);
    for k=1:3
        for l=1:3
            matB(k,l)=b(k)*b(l);
            matC(k,l)=c(k)*c(l);
            localH(k,l,e)=(matB(k,l)+matC(k,l))/(4*D);
        end
    end
    clear ('Xnode','Ynode','node','b','c','matB','matC','D')
end