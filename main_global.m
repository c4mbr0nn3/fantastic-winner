% Topology & coordinates file load
prompt='File Topologia: ';
A=input(prompt,'s');
topol=csvread(A);

prompt2='File Coordinate: ';
B=input(prompt2,'s');
coord=csvread(B);

% Read N, NT, TRIANG, NCOORD
N=coord(1,1);
NT=topol(1,1);
TRIANG=topol(2:NT+1,2:4).';
NCOORD=coord(2:N+1,2:3);
clear('coord','topol','prompt','prompt2','A','B')

% Fix N1=15
N1=15;

% Calculate IA, JA, NTERM
[NTERM,JA,IA]=funtopol(N,NT,N1,TRIANG);

% Calculate TRIJA
[TRIJA]=TRIJA(NT,TRIANG,IA,JA);

% Calculate local Stiffness Matrix
[localH]=localH(TRIANG,NCOORD,NT);

