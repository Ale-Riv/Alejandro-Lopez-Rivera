function tpe=TamPasoEst(Rang, R, bacterias, Var)

[sb,c]=size(bacterias);

   for u=1:sb
       for k=1:Var
          bacterias(u,Var+k)=(Rang(k,1) + (Rang(k,2) - Rang(k,1)))*R;
       end
   end
 tpe=bacterias;
