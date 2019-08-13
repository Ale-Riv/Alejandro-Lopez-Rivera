function bact=evaluacions(bacteria,Var)

[f,c]=size(bacteria);

  eps=0;

    for i=1:f
       x=bacteria(i,1:Var);
       fx=resorte(x);
       [c ceq]=Cresorte(x);
       bacteria(i,(Var+Var+1))=fx;
       [h,y]=size(c);
       svr=0;
          if h>0
             for j=1:h
                 if c(j)>0
                    svr=svr+c(j);
                 end
             end
          end
       [h,y]=size(ceq);
          if h>0
             for j=1:h
                 ceq(j)=abs(ceq(j))-eps;
                 if ceq(j)>0
                    svr=svr+ceq(j);
                 end
             end
          end
      bacteria(i,(Var+Var+2))=svr;
    end

bact=bacteria;
