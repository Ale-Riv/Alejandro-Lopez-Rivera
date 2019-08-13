function tpd=TamPasoDin(Rang, bacterias, Var)

[sb,c]=size(bacterias);

   for u=1:sb
      for k=1:Var
         bacterias(u,Var+k)=(Rang(k,1) + (Rang(k,2) - Rang(k,1)))*rand();
      end
   end

 tpd=bacterias;
           
