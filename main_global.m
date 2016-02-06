%%% NUMERICAL METHODS FOR DIFFERENTIAL EQUATIONS %%%
%%%%%%%%%%%% Zorzi Francesco - 112930 %%%%%%%%%%%%%%
% Solve Poisson PDE with PCG.                      %
% Known term: f= - 2*x^2 - 2*y^2 + 4               %
% Domain: x,y=[-1,1]                               %
% Dirichlet b.c. x=[-1,1],y=-1,1                   %
% Neumann b.c. q=2*(1-y^2),y=[-1,1],x=-1,1         %

clear all

choice=menu('Mesh','0','1','2','3','4');
if choice==1
%     A='mesh0.topol';
%     B='mesh0.coord';
%     C='mesh0.bound';
    A='mesh0.topol';
    B='mesh0.coord';
    C='mesh0.bound';
elseif choice==2
%     A='mesh1.topol';
%     B='mesh1.coord';
%     C='mesh1.bound';
    A='mesh1-renum.topol';
    B='mesh1-renum.coord';
    C='mesh1-renum.bound';
elseif choice==3
%     A='mesh2.topol';
%     B='mesh2.coord';
%     C='mesh2.bound';
    A='mesh2-renum.topol';
    B='mesh2-renum.coord';
    C='mesh2-renum.bound';
elseif choice==4
%     A='mesh3.topol';
%     B='mesh3.coord';
%     C='mesh3.bound';
    A='mesh3-renum.topol';
    B='mesh3-renum.coord';
    C='mesh3-renum.bound';
elseif choice==5
%     A='mesh4.topol';
%     B='mesh4.coord';
%     C='mesh4.bound';
    A='mesh4-renum.topol';
    B='mesh4-renum.coord';
    C='mesh4-renum.bound';
end

choice2=menu('CR Method iteration','0','1','10','20','40','60');
if choice2==1
    maxiter=0;
elseif choice2==2
    maxiter=1;
elseif choice2==3
    maxiter=10;
elseif choice2==4
    maxiter=20;
elseif choice2==5
    maxiter=40;
elseif choice2==6
    maxiter=60;
end

choice3=menu('Preconditioner','Jacobi','Incomplete Cholesky');

% Topology, coordinates and b.cond. file load
topol=csvread(A);
coord=csvread(B);
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
clear('coord','topol','bound','A','B','C')

% Fix N1=15 (max node contacs)
N1=15; 

% Calculate Stiffness Matrix IA, JA, NTERM
[NTERM,JA,IA]=funtopol(N,NT,N1,TRIANG);

% Calculate Stiffness Matrix TRIJA
[TRIJA]=TRIJA(NT,TRIANG,IA,JA);

% Calculate Stiffness Matrix SYSMAT
[SYSMAT,areanod]=stiffmat(TRIANG,NCOORD,NT,NTERM,TRIJA,N);

% Calculate known term f with b.cond.
[b,SYSMAT]=bcond(N,NCOORD,areanod,NDnum,NNnum,NNbound,NDbound,IA,SYSMAT);

if choice3==2
    % Incomplete Cholesky factorization LL^T (Kershaw)
    PREC=kersh(N,NTERM,IA,JA,SYSMAT);
    % Calculate initial point with CR Method
    b=-b;
    [h]=CR(b,N,IA,SYSMAT,maxiter,JA,PREC,NTERM);
elseif choice3==1
    h=zeros(1,N);
    PREC=1;
end

% Solve with PCG (Hu-b=0)
[u,taureal,iter,Rres]=PCG(h,b,SYSMAT,IA,JA,N,NTERM,PREC,choice3);

% Analitical solution, Euclidean error, punctual error
[error,ansol,perror]=euerror(N,NCOORD,areanod,u);