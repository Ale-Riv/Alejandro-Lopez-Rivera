function bact=eliminacionD(Var, Rang)

lim=Rang;
acum = 0.0;
bacterias=zeros(1,Var+Var+2);

angulos=zeros(1,Var);

            for i = 1:Var
                alea=randn();
                while(alea<-0.25 || alea > 0.15)
                    alea=randn();
                end
                angulos(i) = alea;
                acum = acum +(angulos(i)^2.0);
            end

        raiz = sqrt(acum);

     for j=1:Var
           bacterias(1,j)=lim(j,1) + (lim(j,2) - lim(j,1)) * rand(1);
           bacterias(1,Var+j)=angulos(j) / raiz;
     end

bacterias=evaluacions(bacterias,Var);
bact=bacterias;
