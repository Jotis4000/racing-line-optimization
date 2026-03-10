function distBounds = getDistBounds(pos,track)

    % Gets the distance between the current position and the boundaries of
    % the track for each point

    left = pos-track.l
    right = pos-track.r
    distBounds = [sqrt(sum(left.^2,2)) sqrt(sum(right.^2,2))]

end