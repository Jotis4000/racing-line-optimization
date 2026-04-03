function len = getLineLength(pos)

    % Gets the distance of the current racing line

    len = 0;
    % disp(length(pos(:,1)))
    for i = 1:length(pos(:,1))

        if i==length(pos(:,1))
            len = len+norm(pos(1,:)-pos(i,:));
        else
            len = len+norm(pos(i+1,:)-pos(i,:));
        end
        
    end

end