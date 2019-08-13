function [bact]=quimiotaxis(bacteria,Nc,B,Var,Rang)

[f,c]=size(bacteria);
angulos=zeros(Var);
lim=Rang;
mitad=floor(Nc/2);

    for b=1:f
      bandera = 1;
      bacteria(b,Var+Var+3)=0;

        for c=1:Nc   %inicia ciclo quimiotáxico
           newbacteria=zeros(1,Var+Var+2);

             if (bandera == 1)
               angulos=angulo(Var);
             end

           vastago=size(1,Var+Var+2);
           v1=randi(f,1);
           v2=randi(f,1);

             while(v1==v2)
               v2=randi(f,1);
             end

           v3=randi(f,1);

             while(v3==v2 || v3==v1)
                v3=randi(f,1);
             end

           bact1=bacteria(v1,:);
           bact2=bacteria(v2,:);
           bact3=bacteria(v3,:);
           %***************************++

             for m=1:Var
               vastago(1,m)=bact1(1,m)+(B-1)*(bact2(1,m)-bact3(1,m));
             end

             for k=1:Var

                if (c > mitad || c < mitad)
                    if (rem(c,2)==0)
                       newbacteria(1,k)=vastago(1,k);
                       newbacteria(1,k+Var)=bacteria(b,k+Var);
                        else
                       newbacteria(1,k)=bacteria(b,k)+bacteria(b,k+Var)*angulos(k);
                       newbacteria(1,k+Var)=bacteria(b,k+Var);
                    end

                    else
                       %Agrupamiento
                       newbacteria(1,k)=bacteria(b,k)+B*(bacteria(1,k)-bacteria(b,k));
                       newbacteria(1,k+Var)=bacteria(1,k+Var);
                end

                if (newbacteria(1,k) < lim(k,1))
                     newbacteria(1,k) = ((lim(k,1) * 2.0) - newbacteria(1,k));
                    elseif (newbacteria(1,k) > lim(k,2))
                   newbacteria(1,k) = ((lim(k,2) * 2.0) - newbacteria(1,k));
                end

                if (newbacteria(1,k) < lim(k,1))
                   newbacteria(1,k) = lim(k,1) + (lim(k,2) - lim(k,1)) * rand(1);
                    elseif (newbacteria(1,k) > lim(k,2))
                   newbacteria(1,k) = lim(k,1) + (lim(k,2) - lim(k,1)) * rand(1);
                end
            end

           newbacteria=evaluacions(newbacteria,Var);

            if (newbacteria(1,Var+Var+2) == 0 && bacteria(b,Var+Var+2)== 0)
                   if (newbacteria(1,Var+Var+1) < bacteria(b,Var+Var+1))
                      bandera=0;
                      bacteria(b,1:Var+Var+2)=newbacteria(1,:);
                      bacteria(b,Var+Var+3)=bacteria(b,Var+Var+3)+1;
                        else
                       bandera=1;
                   end

              elseif (newbacteria(1,Var+Var+2) > 0 && bacteria(b,Var+Var+2)> 0)
                   if (newbacteria(1,Var+Var+2) < bacteria(b,Var+Var+2))
                      bandera=0;
                      bacteria(b,1:Var+Var+2)=newbacteria(1,:);
                      bacteria(b,Var+Var+3)=bacteria(b,Var+Var+3)+1;
                         else
                       bandera=1;
                   end

              elseif (newbacteria(1,Var+Var+2) ==0 &&  bacteria(b,Var+Var+ 2)> 0)
                    bandera = 0;
                    bacteria(b,1:Var+Var+2)=newbacteria(1,:);
                    bacteria(b,Var+Var+3)=bacteria(b,Var+Var+3)+1;
              else
                  bandera=1;
           end
       end
      bacteria=ordenar(bacteria,Var);
    end
bact=bacteria;
