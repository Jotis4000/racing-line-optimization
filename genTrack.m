function track = genTrack()

    %% Parameters
    addpath("tracks\")

    FileName = "Zandvoort_centerline.csv";
    datat = xlsread(FileName)*10; % x_m, y_m, w_right, w_left

    track.x_m = datat(:,1)
    track.y_m = datat(:,2)

    % Generate track bounds - iterate through data set, find direction,
    % add point perpendicular
    % There is currently a redundant point at the origin. No clue how to
    % fix it or what to do

    for i = 1:length(track.x_m)
        
        if i==length(track.x_m)
            vec = [track.x_m(1)-track.x_m(i) track.y_m(1)-track.y_m(i)]
        else
            vec = [track.x_m(i+1)-track.x_m(i) track.y_m(i+1)-track.y_m(i)]
        end

        mag = norm(vec)

        disp(mag)
        % phi = 

    end

    plot(track.x_m,track.y_m)
    % disp(datat(:,1))

end