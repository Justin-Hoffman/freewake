function plotState(surfs)
        alpha = 0.3;
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
            hSurface = surfl(surfs(n).xSurface,-surfs(n).ySurface,-surfs(n).zSurface);
            material dull
            set(hSurface,'FaceColor',[.5 .5 .5],'FaceAlpha',0.9);
            hSurface = surfl(surfs(n).xWake,-surfs(n).yWake,-surfs(n).zWake);
            set(hSurface,'FaceColor',cWake{n},'FaceAlpha',alpha);
            material dull
            camlight right
            camlight(30,80)
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
            plot(surfs(n).yCp(:,1),-surfs(n).zSpanwiseForce)
            xlabel('Y');
            ylabel('Sectional Loading');
            grid on
        end

        figure(103);
        subplot 311
        hold off
        plot(surfs(1).T, -surfs(1).CFZ/2); %Divide by 2 for rotor CT
        grid on
        xlabel('Time'); ylabel('CT');
        subplot 312 
        hold off
        plot(surfs(1).T, -surfs(1).CMZ/2); %Divide by 2 for rotor CT
        grid on
        xlabel('Time'); ylabel('CQ');
        subplot 313 
        hold off
        plot(surfs(1).T, (-surfs(1).CFZ/2).^(3/2)./(sqrt(2)*-surfs(1).CMZ/2)); %Divide by 2 for rotor CT
        grid on
        xlabel('Time'); ylabel('FM');
        drawnow;

        
    end