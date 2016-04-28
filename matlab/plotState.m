function plotState(surfs)
        
        nSurfs = length(surfs);
        cWake = {'r','b','g'};
        %ColorSet = varycolor(Nb*(N_span+1));
        figure(101); clf; hold all;
        %set(gca, 'ColorOrder', ColorSet);
        for n = 1:nSurfs
           sSurf = size(surfs(n).xSurface);
           sWake = size(surfs(n).xWake);
           N_span = sSurf(1);
           N_chord = sSurf(2);
           N_wake = sWake(2);
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
%            for i = 1:N_span
%               plot3(surfs(n).xWake(i,:),...
%                     -surfs(n).yWake(i,:),...
%                     -surfs(n).zWake(i,:), cWake{n} )
%            end
%            for j = 1:N_wake
%               plot3(surfs(n).xWake(:,j),...
%                     -surfs(n).yWake(:,j),...
%                     -surfs(n).zWake(:,j), cWake{n} )
%            end
        end
        axis equal;
        view([0,1,.1]);
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