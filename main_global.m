% Define known term f (our problem: - 2*x^2 - 2*y^2 + 4)
syms x y
prompt0='Equazione f: ';
f(x,y)=input(prompt0);

% Topology & coordinates file load
prompt='File Topologia: ';
A=input(prompt,'s');
topol=csvread(A);

prompt2='File Coordinate: ';
B=input(prompt2,'s');
coord=csvread(B);

% Boundary condition nodes
prompt3='File B.C.: ';
C=input(prompt3,'s');
bound=csvread(C);

% Store N, NT, TRIANG, NCOORD
N=coord(1,1);
NT=topol(1,1);
TRIANG=topol(2:NT+1,2:4).';
NCOORD=coord(2:N+1,2:3);
NDbound=bound(1,:);
NNbound=bound(2,:);
NDbound=NDbound(NDbound~=0);
NNbound=NNbound(NNbound~=0);
clear('coord','topol','prompt','prompt2','A','B','C','bound')

% Fix N1=15
N1=15;

% Calculate Stiffness Matrix IA, JA, NTERM
[NTERM,JA,IA]=funtopol(N,NT,N1,TRIANG);

% Calculate Stiffness Matrix TRIJA
[TRIJA]=TRIJA(NT,TRIANG,IA,JA);

% Calculate Stiffness Matrix SYSMAT
[SYSMAT,areanod]=stiffmat(TRIANG,NCOORD,NT,NTERM,TRIJA,N);

% Calculate known term f
b=zeros(1,N);
for k=1:N
    Xnode=NCOORD(k,1);
    Ynode=NCOORD(k,2);
    b(k)=f(Xnode,Ynode)*areanod(k);
    clear('Xnode','Ynode')
end
