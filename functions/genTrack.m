function track = genTrack(plt)

    %% Parameters
    addpath("tracks\")

    FileName = "Budapest.csv";
    datat = xlsread(FileName); % x_m, y_m, w_right, w_left

    track.m = [datat(:,1) datat(:,2)]; % X, Y
    track.w = [datat(:,4) datat(:,3)]; % Left, Right

    % Generate track bounds - iterate through data set, find direction,
    % add point perpendicular

    track.vec = [circshift(track.m(:,1), [-1 0])-track.m(:,1) circshift(track.m(:,2), [-1 0])-track.m(:,2)];
    track.vecmag = vecnorm(track.vec,2,2);
    track.l = track.m+[-track.vec(:,2), track.vec(:,1)]./track.vecmag.*track.w(:,1);
    track.r = track.m+[track.vec(:,2), -track.vec(:,1)]./track.vecmag.*track.w(:,2);

    if plt
        plot(track.m(:,1),track.m(:,2),track.l(:,1),track.l(:,2),track.r(:,1),track.r(:,2))
        legend("m","l","r")
    end

    track.mu = 1.4; % Ground friction coefficient

end