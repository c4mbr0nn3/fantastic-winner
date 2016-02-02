function [NTERM,JA,IA]=funtopol(N,NT,N1,TRIANG)
 NTERM=N1*N;
 JA=zeros(1,NTERM);
 JA(1)=1;
 for K=1:(NTERM-1)
      if (fix(K/N1)*N1== K )
         JA(K+1)=K/N1+1;
      end
end
%
%  analizza tutti i triangoli ( NT )
%
for  K=1:NT
     for I=1:3
          I1(I)=TRIANG(I,K);
     end
%
%  ordina i nodi dell'elemento in senso crescente 
%
     for J=1:3
         for I=1:2
               if ( I1(I) > I1(I+1) )
                  aux=I1(I);
		  I1(I)=I1(I+1);
		  I1(I+1)=aux;
	       end 
         end
     end 
%
%  alla fine del ciclo in I1 abbiamo la terna ordinata 
%
%  
%  genera il vettore JA
%
   for I=1:2
       J=I+1;
       for  L=J:3
           M=N1*(I1(I)-1)+L-J+2;
           MCONTR=N1*I1(I);
           continua=1;
% ciclo while # 1	      
           while ( (((I1(L)-JA(M))<0) | ((I1(L)-JA(M))> 0)) & (continua==1)  )
              if ((I1(L)-JA(M))>0)
                 if (JA(M)==0 )
                     JA(M)=I1(L);
                     continua=0;
                  else
                      M=M+1;
                      if ((M-MCONTR)>=0)
                         ivar3=MCONTR/N1;
                         K, M, MCONTR
                         error('MATLAB:topologia', 'errore ')
                       else
                          continua=1;
                       end 
                   end 
                elseif ((I1(L)-JA(M))<0)
                      MM=M;
                      MM=MM+1;
                      if (M-MCONTR>=0)
                          ivar3=MCONTR/N1;
                          K
                          error('MATLAB:topologia', 'errore ')
                       else
                       while (JA(MM)>0 | JA(MM)<0)
                           MM=MM+1;
                           if ( M-MCONTR>=0)
                              ivar3=MCONTR/N1;
                              K
                              error('MATLAB:topologia', 'errore ')
			   end
                       end 
                       while (JA(MM)==0 )
                           JA(MM)=JA(MM-1);
                           MM=MM-1;
                           while (MM>M )
                              JA(MM)=JA(MM-1);
                              MM=MM-1;
                          end 
                        end 
                 end   %fine ciclo if
                 JA(M)=I1(L);
             end  % fine ciclo if
        end %fine ciclo while # 1
      end   % fine ciclo for L=J:3
   end  % fine ciclo for i=1:3
end % fine ciclo  for K=1:NT
                         
%
%  costruisce il vettore IA
%
IA(1)=1;
M=1;
J=1;
K=1;
 for K=1:N1:NTERM
      for I=1:N1
          if (JA(K+I-1) ~= 0)
            M=M+1;
          end 
      end
      J=J+1;
      IA(J)=M;
 end
%
%  compatta il vettore JA eliminando gli zeri
%
 M=0;
 for K=1:NTERM
     if (JA(K)~=0)
          M=M+1;
          JA(M)=JA(K);
     end 
end
 NTERM=M;
 JA=JA(1:NTERM);

