function h=pcolor2(x,y,z)
% pcolor2(x,y,z): pcolor front end that plots everything and assumes
% x and y are the mid-point locations of the grid boxes.
%
% Example:
%   [x,y]=meshgrid(3:7,1:4); z=rand(size(x));
%   subplot(2,1,1), pcolor(x,y,z), shading flat
%   subplot(2,1,2), pcolor2(x,y,z)
%
% Ian Eisenman, 2006

if nargin<3 % x and y not given
    % remove singleton dimensions
    x=squeeze(x);
    z=x;
    [x,y]=meshgrid(1:size(z,2),1:size(z,1));
end

if min(size(x))==1 % x,y are vectors, not matrices
    [x,y]=meshgrid(x,y);
end

% remove singleton dimensions
z=squeeze(z);

% make sure x is rows and y is columns
if abs(x(2,1)-x(1,1)) > abs(x(1,2)-x(1,1)) % x increases downward
    x=x';
    y=y';
    z=z';
end

dx=diff(x(1,:));
dy=diff(y(:,1));
x=[x; x(end,:)]; x=[x x(:,end)+dx(end)]-repmat(dx([1 1:end end]),[size(x,1),1])/2;
y=[y; y(end,:)+dy(end)]; y=[y y(:,end)]-repmat(dy([1 1:end end]),[1,size(y,2)+1])/2; 

% duplicate last row and column so that pcolor will show them
z=[z; z(end,:)]; z=[z z(:,end)]; 

h0=pcolor(x,y,z);

shading flat;

if(nargout > 0)
    h = h0;
end
