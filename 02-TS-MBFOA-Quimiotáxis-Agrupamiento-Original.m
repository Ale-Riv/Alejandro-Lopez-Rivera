function [bact]=quimiotaxis(bacteria,prob,nc)
[f,c]=size(bacteria);
var=numVariables(prob);
angulos=zeros(var);
lim=rangoVariables(prob);
mitad=floor(nc/2);
factorE=1.75;
        for b=1:f
            bandera = 1;
            bacteria(b,var+var+3)=0;
            for c=1:nc   %inicia ciclo quimiotáctico
                newbacteria=zeros(1,var+var+2);
                if (bandera == 1)
                    angulos=angulo(var);
                end

                vastago=size(1,var+var+2);
                v1=randi(f,1);
                v2=randi(f,1);
                while(v1==v2)
                     v2=randi(f,1);
                end
                v3=randi(f,1);
                while(v3==v2 || v3==v1)
                     v3=randi(f,1);
                end

                bact1=bacteria(b,:);

         %else
          % bact1=bacterias(v1,:);
         %end

                bact2=bacteria(v2,:);
                bact3=bacteria(v3,:);
                %***************************++

          for m=1:var
              vastago(1,m)=bact1(1,m)+(factorE-1)*(bact2(1,m)-bact3(1,m));
          end
         % vastago=bact1+angulos(1)*(bact2-bact3);

           %vastago2=evaluacions(vastago,prob);

           %if(vastago2(var+var+2)==0 && (vastago2(var+var+1)<bacteria(b,var+var+1)))
            %   bacteria(b,1:var)=vastago2(1:var);
             %  bacteria(b,var+var+1:var+var+2)=vastago2(var+var+1:var+var+2);

           %end





               for k=1:var
                        %if(c==1000000000000000)
                      % if (rem(gene,100)==0 && c==fsw &&( b==floor(f/2)))

                     %  factorE=rand(1)+1;
                      %  newbacteria(1,k) = bacteria(b,k)+factorE*(bacteria(1,k)-bacteria(b,k));
                       %newbacteria(1,k)=bacteria(1,k) + 0.7 *(bact2(1,k)-bact3(1,k));
                      %  newbacteria(1,k)=bact1(1,k) + 0.7 *(bact2(1,k)-bact3(1,k));
                       % else
                        %% if(rem(c,2)==0)
                            %newbacteria(1,k) = bacteria(b,k)+bacteria(b,k+var)*angulos(k);
                            % newbacteria(1,k) = bacteria(b,k)+bacteria(b,k+var);
                          %newbacteria(1,k) = bacteria(b,k)+angulos(k);
                  %%newbacteria(1,k)=bact1(1,k) + 0.7 *(bact2(1,k)-bact3(1,k));

                  %%   else
                              % newbacteria(1,k) = bacteria(b,k)+(angulos(k)*bacteria(b,k+var));
                          %newbacteria(1,k) = bacteria(b,k)+(bacteria(b,k+var)*angulos(k));
                          %version1
                          %newbacteria(1,k)=bacteria(b,k)+ angulos(k)* (bacteria(b,k+var) *(bact2(1,k)-bact3(1,k)));
                         %version2
                  %%newbacteria(1,k)=bacteria(b,k)+ (bacteria(b,k+var) *(bact2(1,k)-bact3(1,k)));
                          % newbacteria(1,k)=bacteria(b,k) + bacteria(b,k)*bacteria(b,k+var);
                         %newbacteria(1,k) = bacteria(b,k)+(0.001*angulos(k));
                  %%  end
                       % end
                     %newbacteria(1,k+var)=angulos(k)*0.5;
                      %newbacteria(1,k+var)=bacteria(b,k+var)*angulos(k);
                      if (c > mitad || c < mitad)
                          if (rem(c,2)==0)
                         newbacteria(1,k)=vastago(1,k);
                          newbacteria(1,k+var)=bacteria(b,k+var);
                          else
                            newbacteria(1,k)=bacteria(b,k)+bacteria(b,k+var)*angulos(k);
                            newbacteria(1,k+var)=bacteria(b,k+var);
                          end


                      else
                          newbacteria(1,k)=bacteria(b,k)+factorE*(bacteria(1,k)-bacteria(b,k));
                           newbacteria(1,k+var)=bacteria(1,k+var);

                      end

                        if (newbacteria(1,k) < lim(k,1))
                            newbacteria(1,k) = ((lim(k,1) * 2.0) - newbacteria(1,k));
                        else if (newbacteria(1,k) > lim(k,2))
                            newbacteria(1,k) = ((lim(k,2) * 2.0) - newbacteria(1,k));
                            end
                        end
                        if (newbacteria(1,k) < lim(k,1))
                            newbacteria(1,k) = lim(k,1) + (lim(k,2) - lim(k,1)) * rand(1);
                        else if (newbacteria(1,k) > lim(k,2))
                            newbacteria(1,k) = lim(k,1) + (lim(k,2) - lim(k,1)) * rand(1);
                            end
                        end
                end

               newbacteria=evaluacions(newbacteria,prob);

                if (newbacteria(1,var+var+2) == 0 && bacteria(b,var+var+2)== 0)
                    if (newbacteria(1,var+var+1) < bacteria(b,var+var+1))
                        bandera=0;
                            bacteria(b,1:var+var+2)=newbacteria(1,:);
                            bacteria(b,var+var+3)=bacteria(b,var+var+3)+1;
                    else
                        bandera=1;
                    end
                elseif (newbacteria(1,var+var+2) > 0 && bacteria(b,var+var+2)> 0)
                    if (newbacteria(1,var+var+2) < bacteria(b,var+var+2))
                        bandera=0;
                        bacteria(b,1:var+var+2)=newbacteria(1,:);
                         bacteria(b,var+var+3)=bacteria(b,var+var+3)+1;
                    else
                        bandera=1;
                    end
                elseif (newbacteria(1,var+var+2) ==0 &&  bacteria(b,var+var+ 2)> 0)
                    bandera = 0;
                   bacteria(b,1:var+var+2)=newbacteria(1,:);
                    bacteria(b,var+var+3)=bacteria(b,var+var+3)+1;
                else
                    bandera=1;
                end
            end
              bacteria=ordenar(bacteria,prob);
        end
bact=bacteria;
