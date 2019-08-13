function bact=eliminacionD( prob)
    var=numVariables(prob);
    lim=rangoVariables(prob);
    acum = 0.0;
    bacterias=zeros(1,var+var+2);
angulos=zeros(1,var);

for i = 1:var
    alea=randn();
    while(alea<-0.25 || alea > 0.15)
        alea=randn();
    end
    angulos(i) = alea;
            acum = acum +(angulos(i)^2.0);
end
    raiz = sqrt(acum);

     for j=1:var
           bacterias(1,j)=lim(j,1) + (lim(j,2) - lim(j,1)) * rand(1);
           bacterias(1,var+j)=angulos(j) / raiz;
     end

bacterias=evaluacions(bacterias,prob);
  bact=bacterias;
