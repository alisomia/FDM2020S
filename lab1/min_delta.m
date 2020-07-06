function s = min_delta(x,y)
global delta
    s = (x + y - sqrt((x-y).^2+delta))/2;
    