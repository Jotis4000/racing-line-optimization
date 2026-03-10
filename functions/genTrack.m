function track = genTrack()

    %% Parameters
    addpath("tracks\")

    FileName = "Spa.csv";
    datat = xlsread(FileName)*10; % x_m, y_m, w_right, w_left

    track.m = [datat(:,1) datat(:,2)] % X, Y
    track.w = [datat(:,4) datat(:,3)] % Left, Right

    track.l = zeros(length(track.m(:,1)),2)
    track.r = zeros(length(track.m(:,1)),2)

    % Generate track bounds - iterate through data set, find direction,
    % add point perpendicular

    for i = 1:length(track.m(:,1))

        if i==length(track.m(:,1))
            vec = [track.m(1,1)-track.m(i,1) track.m(1,2)-track.m(i,2)]
        else
            vec = [track.m(i+1,1)-track.m(i,1) track.m(i+1,2)-track.m(i,2)]
        end

        mag = norm(vec)
        Rright = [0 -1; 1 0] % I thought these were backwards? lol?
        Rleft = [0 1; -1 0]

        vecleft = vec*Rleft/mag
        vecright = vec*Rright/mag

        % disp(vecleft)

        track.l(i,1) = track.m(i,1)+vecleft(1)*track.w(i,1)
        track.l(i,2) = track.m(i,2)+vecleft(2)*track.w(i,1)

        track.r(i,1) = track.m(i,1)+vecright(1)*track.w(i,2)
        track.r(i,2) = track.m(i,2)+vecright(2)*track.w(i,2)

        % disp(mag)
        % phi = 

    end

    % plot(track.x_m,track.y_m,track.x_l,track.y_l,track.x_r,track.y_r)
    plot(track.m(:,1),track.m(:,2),track.l(:,1),track.l(:,2),track.r(:,1),track.r(:,2))
    legend("m","l","r")
    % disp(datat(:,1))

end