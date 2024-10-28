function y = findfirstpoint(y1,y2)

    y = 2*y1 - y2;

    if y2 < y1

        if y < y1
            y = y1 + abs(y1 - y2);
        end

    else

        if y > y1
            y = y1 - abs(y1 - y2);
        end

    end