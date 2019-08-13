function bact=ordenar(bacteria, Var)

[f,c]=size(bacteria);
aux= zeros(1,Var+Var+2);

        for i=1:f
            for j=1:(f-1)
                if (bacteria(j,(Var+Var+2))==0 && bacteria(j+1,(Var+ Var+ 2)) == 0)
                    if (bacteria(j,(Var+Var+1)) > bacteria(j + 1,(Var+ Var+1)))
                        for k=1:c
                            aux(k) = bacteria(j,k);
                            bacteria(j,k) = bacteria(j + 1,k);
                            bacteria(j + 1,k) = aux(k);
                        end
                    end
                end
                if (bacteria(j,(Var+Var+2))>0 && bacteria(j+1,(Var+ Var+ 2))>0)
                    if (bacteria(j,(Var+Var+2)) > bacteria(j + 1,(Var+ Var+2)))
                        for k=1:c
                            aux(k) = bacteria(j,k);
                            bacteria(j,k) = bacteria(j + 1,k);
                            bacteria(j + 1,k) = aux(k);
                        end
                    end
                end
                if (bacteria(j,(Var+Var+2))> 0 && bacteria(j+1,(Var+ Var+ 2)) == 0)
                        for k=1:c
                            aux(k) = bacteria(j,k);
                            bacteria(j,k) = bacteria(j + 1,k);
                            bacteria(j + 1,k) = aux(k);
                        end
                end
            end
       end

bact=bacteria;
