function bact=eliminacionD(num, prob)
    var=numVariables(prob);
    lim=rangoVariables(prob);

    bacterias=zeros(num,var+var+2);


for i=1:num
         acum = 0.0;
        angulos=zeros(1,var);
 for y = 1:var
    alea=randn();
    while(alea<-1 || alea > 1)
        alea=randn();
    end
             angulos(y) = alea;
            acum = acum +(angulos(y)^2.0);
 end

     for j=1:var
          bacterias(i,j)=lim(j,1) + (lim(j,2) - lim(j,1)) * rand(1);
          bacterias(i,var+j)=lim(j,1) + (lim(j,2) - lim(j,1)) * rand(1);
     end


end

bacterias=evaluacions(bacterias,prob);
%bacterias=normalizacion(bacterias,tipoN, prob);
  bact=bacterias;
  
