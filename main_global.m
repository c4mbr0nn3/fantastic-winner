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

% Store N, NT, TRIANG, NCOORD, NDnum, NNnum, NDbound, NNbound
N=coord(1,1); 
NT=topol(1,1); 
TRIANG=topol(2:NT+1,2:4).'; 
NCOORD=coord(2:N+1,2:3); 
NDnum=bound(1,1);
NNnum=bound(2+ceil(NDnum/10),1);
NDbound=bound(2:1+ceil(NDnum/10),1:10).'; 
NDbound=NDbound(NDbound~=0).';
NNbound=bound(3+ceil(NDnum/10):2+ceil(NDnum/10)+ceil(NNnum/10),1:10).'; 
NNbound=NNbound(NNbound~=0).'; 
clear('coord','topol','bound','prompt','prompt2','prompt3','prompt0','A','B','C')

% Fix N1=15
N1=15;

% Calculate Stiffness Matrix IA, JA, NTERM
[NTERM,JA,IA]=funtopol(N,NT,N1,TRIANG);

% Calculate Stiffness Matrix TRIJA
[TRIJA]=TRIJA(NT,TRIANG,IA,JA);

% Calculate Stiffness Matrix SYSMAT
[SYSMAT,areanod]=stiffmat(TRIANG,NCOORD,NT,NTERM,TRIJA,N);

% Calculate known term f with b.cond.
[b,SYSMAT]=bcond(N,NCOORD,areanod,NDnum,NNnum,NNbound,NDbound,IA,SYSMAT);

% Incomplete Cholesky factorization LL^T (Kershaw)
PREC=kersh(N,NTERM,IA,JA,SYSMAT);

% calculate initial point (w=x0,v=b)
v=-b;
[w]=avind(N,PREC,IA,JA,NTERM,v);

% Solve with PCG
h=w;
b=-b;
[x,tau,taureal,iter]=PCG(h,b,SYSMAT,IA,JA,N,NTERM,PREC);

% Plot some shit
TR = triangulation(TRIANG.',[NCOORD(:,1),NCOORD(:,2),x.']);
trimesh(TR)
% semilogy(1:1:iter,tau,'*')
