%% 3D Animation
global index_view;
global old_position;

switch video
    case 'play'
        %set-up
        droneFigure = figure('name','Plot Trajectory','units','pixels','position',[0,0,1280,720]);
        ha = axes('Parent',droneFigure);
        set(droneFigure,'menubar','figure','renderer','opengl');
        set(ha,'Visible','On','Box','On','XGrid', 'on','YGrid', 'on','ZGrid',...
            'on','projection','perspective');
        hold on;
        cameratoolbar('show');
        axis vis3d;
        view(3);
        zoom(0.9);
        index_view = 0;
        old_position = [0 0 0];
        cameratoolbar('ResetCameraAndSceneLight');
        
        %play
        for i = 1:length(sim_out.time)
            draw_quad([sim_out.pos_ned(i,:), sim_out.attitude(i,:), sim_out.motors(i,:)], drone.arm_length, drone.propeller_radius);
            pause(1/init.fps);
        end
        drawnow;
        
    case 'playAndRecord'
        %set-up
        droneFigure = figure('name','Plot Trajectory','units','points','position',[0,0,640,480]);
        ha = axes('Parent',droneFigure);
        set(droneFigure,'menubar','figure','renderer','opengl');
        set(ha,'Visible','On','Box','On','XGrid', 'on','YGrid', 'on','ZGrid',...
            'on','projection','perspective');
        hold on;
        cameratoolbar('show');
        axis vis3d;
        view(3);
        zoom(0.9);
        index_view = 0;
        old_position = [0 0 0];
        cameratoolbar('ResetCameraAndSceneLight');
        
        quadmovie = VideoWriter(videoName);
        quadmovie.FrameRate = init.fps;
        open(quadmovie);
        
        %play and record
        for i = 1:length(sim_out.time)
            draw_quad([sim_out.pos_ned(i,:), sim_out.attitude(i, :), sim_out.motors(i,:)], drone.arm_length, drone.propeller_radius);
            F = getframe(gcf);
            writeVideo(quadmovie,F);
        end
        drawnow;close
        close(quadmovie);
        clear F
        
    otherwise
        %nothing to do
end

clear index_view old_position video videoName i ha