function s = min_delta_grad(x,y)
global delta
s = (1 - (x-y)./sqrt((x-y).^2+delta))/2;