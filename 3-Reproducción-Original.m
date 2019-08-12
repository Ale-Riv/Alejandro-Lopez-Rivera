function bact=reproduccion(bacteria, repro)
[f,c]=size(bacteria);
        for y=(f-repro)+1:f
            for h=1:c
                bacteria(y,h)=bacteria((f-y)+1,h);
            end
        end
bact=bacteria;
