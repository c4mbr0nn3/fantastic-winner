% Global Stiffness Matrix function
function [SYSMAT,areanod]=stiffmat(TRIANG,NCOORD,NT,NTERM,TRIJA,N)
SYSMAT=zeros(1,NTERM); % Initialize SYSMAT
areanod=zeros(1,N); % Initialize areanod
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
    b(3)=Ynode(1)-Ynode(2);
    c(1)=Xnode(3)-Xnode(2);
    c(2)=Xnode(1)-Xnode(3);
    c(3)=Xnode(2)-Xnode(1);
    D=det([ones(1,3);Xnode;Ynode].')/2;
    
    % Calculate areanod for further application
    for k=1:3
        l=node(k);
        areanod(l)=areanod(l)+D/3;
    end
    
    % Initialize useful matrix
    matB=zeros(3,3);
    matC=zeros(3,3);
    localH=zeros(3,3);
    
    % Calculate Stiffness Matrix
    for k=1:3
        for l=1:3
            matB(k,l)=b(k)*b(l);
            matC(k,l)=c(k)*c(l);
            localH(k,l)=(matB(k,l)+matC(k,l))/(4*D); % si può sfruttare simmetria della matrice H locale?
            if TRIJA(k,l,e)~=0
                ind=TRIJA(k,l,e);
                SYSMAT(ind)=SYSMAT(ind)+localH(k,l);
            end
        end
    end
    clear ('Xnode','Ynode','node','b','c','matB','matC','D','localH')
end