function plotState(surfs)
        
        nSurfs = length(surfs);
        cWake = {'r','b','g'};
        %ColorSet = varycolor(Nb*(N_span+1));
        figure(101); clf; hold all;
        %set(gca, 'ColorOrder', ColorSet);
        for n = 1:nSurfs
           sSurf = size(surfs(n).xSurface);
           sWake = size(surfs(n).xWake);
           sFilaments = size(surfs(n).xTipFilament);
           N_span = sSurf(1);
           N_wake = sWake(2);
           N_fil = sFilaments(2);
%            for i = 1:N_span
%               plot3(surfs(n).xSurface(i,:),...
%                     -surfs(n).ySurface(i,:),...
%                     -surfs(n).zSurface(i,:),'-k')
%            end
%            for j = 1:N_chord
%               plot3(surfs(n).xSurface(:,j),...
%                     -surfs(n).ySurface(:,j),...
%                     -surfs(n).zSurface(:,j),'-k')
%            end
            hSurface = surf(surfs(n).xSurface,-surfs(n).ySurface,-surfs(n).zSurface);
            set(hSurface,'FaceColor',[.5 .5 .5],'FaceAlpha',0.9);
            hSurface = surf(surfs(n).xWake,-surfs(n).yWake,-surfs(n).zWake);
            set(hSurface,'FaceColor',cWake{n},'FaceAlpha',0.2);
           for i = 1:2
              plot3(surfs(n).xTipFilament(i,:),...
                    -surfs(n).yTipFilament(i,:),...
                    -surfs(n).zTipFilament(i,:), cWake{n} )
           end
        end
        axis equal;
        view([1,0,.1]);
        drawnow;
            
        figure(102); clf; hold all;
        for n = 1:nSurfs
            plot(surfs(n).yCp,-surfs(n).zSpanwiseForce)
        end

        figure(103); 
        hold off
        plot(surfs(1).T, -surfs(1).CT/2); %Divide by 2 for rotor CT
        xlabel('Time'); ylabel('CT');
        drawnow;
        
    end