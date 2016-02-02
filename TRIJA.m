function [TRIJA]=TRIJA(NT,TRIANG,IA,JA)
TRIJA=zeros(3,3,NT);
for k=1:NT
    I2=TRIANG(:,k);
    for i=1:3
        ii=I2(i);
        for j=1:3
            jj=I2(j);
            TRIJA(i,j,k)=0;
            if jj>=ii
                ind=IA(ii);
                while (JA(ind)-jj<0 || JA(ind)-jj>0)
                    ind=ind+1;
                end
                TRIJA(i,j,k)=ind;
            end
        end
    end
end