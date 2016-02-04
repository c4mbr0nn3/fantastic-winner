      function PREC=kersh(N,NTERM,IA,JA,SYSMAT)
% traduzione in MATLAB della subroutine scritta in FORTRAN kersh.f
%
%  OSSERVAZIONI UTILI per il passaggio da un codice FORTRAN ad uno MATLAB
% Nel passare dalle istruzioni FORTRAN a quelle MATLAB occorre tener presente:
% + in FORTRAN l'istruzione  ELSE IF si puo' scrivere come ELSE IF o ELSEIF
%   in MATLAB  l'analoga istruzione deve essere scritta tutta attaccata
%    ELSEIF altrimenti si apre un altro ciclo IF
% + in FORTRAN un ciclo del tipo DO k=i,j ... END DO dichiara la variabile
%   k come variabile intera e, appena entra nel ciclo la pone uguale a i
%   k=i. Se j<i il ciclo non da' risultati in uscita ma la variabile k
%   rimane comunque del valore uguale a i.
%   Invece, quando i<j, il ciclo DO va avanti e k viene incrementato  di
%   una unita'. Quando k=j, si "lavora" per l'ultima volta dentro il ciclo
%   perche' poi k viene ancora incrementato di un'unita' e quindi il 
%   suo valore diventa k=j+1 >j e  si esce dal ciclo. Il valore in
%   uscita di k e' dunque j+1
%   In MATLAB, invece, una volta che si esce da un ciclo do k=i:j
%   la variabile k e' come cancellata e quindi non puo' essere utilizzata
%   come variabile all'esterno del ciclo.
%   Da qui si puo' osservare come la variabile k2, nel programma FORTRAN
%   e' usata all'esterno del ciclo do k2=i,j-1 con il valore k2=j
%   (sia se i=j sia per i<j)
%         if(j.ge.i) then
%            prec(ia(ja(j))) = prec(ia(ja(j))) + prec(k2)**2
%        end if
%   Lo stesso discorso non si puo' fare in MATLAB perche' k2 viene cancellato
%   al di fuori del ciclo ma le istruzioni analoghe diventano
%         if  j >= i
%              prec(ia(ja(j))) = prec(ia(ja(j))) + prec(j)^2
%        end
%
%
%
%
%  nequ   : numero di equazioni del sistema (n)
%  nterm  : numero di coefficienti non nulli della matrice (nt)
%  ia     : vettore topologico 
%  ja     : vettore con gli indici di colonna
%  sysmat : vettore per la memorizzazione compatta della matrice del sistema
%
%
PREC=zeros(1,NTERM);

      for kk=1:N-1;
         k = IA(kk);
         a = SYSMAT(k) - PREC(k);
         if a<=0
            display('attenzione')
             kk,k,a
             PREC(IA(kk-1))
            a = (PREC(IA(kk-1)))^2;
         end 
         PREC(k) = sqrt(a);
 
         i = IA(kk) + 1;
         j = IA(kk+1) - 1;


         for k1 = i:j;
            PREC(k1) = (SYSMAT(k1)-PREC(k1))/PREC(k);
         end 

         for k2 = i:j-1;
            j1 = IA(JA(k2));
            PREC(j1) = PREC(j1) + PREC(k2)^2;
            i1 = k2 + 1;
            j1 = j1 + 1;
            while (j1<IA(JA(k2)+1) && i1<=j);
               if(JA(j1)==JA(i1));
                  PREC(j1) = PREC(j1) + PREC(k2)*PREC(i1);
                  i1 = i1 + 1;
                  j1 = j1 + 1;
               elseif (JA(j1)<JA(i1));
                  j1 = j1 + 1;
               elseif (JA(j1)>JA(i1));
                   i1 = i1 + 1;
               end
            end 
         end 

         if  j >= i  ;
               PREC(IA(JA(j))) = PREC(IA(JA(j))) + PREC(j)^2;
         end 

      end 

      k = IA(N);
      a = SYSMAT(k)-PREC(k);
      if(a<=0) 
        display('attenzione')
          N ,a
          PREC(IA(N-1))
         a = (PREC(IA(N-1)))^2;
      end 
      PREC(k) = sqrt(a);

