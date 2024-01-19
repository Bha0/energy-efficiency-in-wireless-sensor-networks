function varargout=drawcircle(center,radius,varargin)
global t DS DE
% This code is draw circle with different center,radius,type and positions
%The input
% center co-ordinated   - [r c];
% radius                -    r
% Boundary              -   [LD UD] LD=Lower limit of Degree, UD upper limit of the degree
% 'quarter'              -   Draw the quater circle
% 'half'                -   Draw half circle
% 'full'                -   Draw full circle
% 'clk'                 -   Draw circle clock wise
% 'aclk'                -   Draw circle Anti clock wise
% '1P'                  -   1st quardrant (+,+)
% '2P'                  -   2nd quardrant (-,+)
% '3P'                  -   3rd quardrant (-,-)
% '4P'                  -   4th quardrant (+,-)

% Example
    % [row col]=drawcircle([1 1],2,'quarter','clk','1P');
    %   Draw the quarter circle with center(1,1), radius-2,
    %   1st quardrant and circle draw clock wise position
    % [row col]=drawcircle([1 1],2,[-90 90]);
    %   Draw circle with spefic Degree range
    % drawcircle([1 1],2,[-90 90]);
    %   Draw only circle
% Created by Sundha R 28.02.2014 
% Based on the basic formula circle with co-ordinates and radius
if length(varargin)== 1
    D   =   varargin{1};
    if D(1) >=  D(2)
        error('The Low limit must be smaller than the upper limit')
        return
    end
    type  =  ' ';
    CF    =  'clk';
elseif length(varargin) ==  2
    D = varargin{1};
    if D(1) >=  D(2)
        error('The Low limit of the Degree must smaller than the uper limit')
        return
    end
    type =  ' ';
    CF   =  varargin{2};
elseif  length(varargin) == 3
    type = varargin{1};
    CF   = varargin{2};
    PO   = varargin{3};
else
    error('Check the input Argument')
end
switch type
    case 'quarter'
        switch PO
            case '1P'
                t = pi:0.0067*pi:(3*pi/2);
            case '2P'
                t = (pi/2):0.0067*pi:pi;
            case '3P'
                t = 0:0.0067*pi:(pi/2);
            case '4P'
                t = (3*pi/2):0.0067*pi:2*pi;
         end
     case 'half'
         switch PO
             case '1P'
                t = pi:0.0067*pi:2*pi;
             case '2P'
                t = (pi/2):0.0067*pi:(3*pi/2);  
             case '3P'
                t = 0:0.0067*pi:pi; 
             case '4P'
                t = (-pi/2):0.0067*pi:(pi/2); 
          end
    case 'full'
        switch PO
        case '1P'
                t = 0:0.0067*pi:2*pi;
          case '2P'
                t = 0:0.0067*pi:2*pi;
          case '3P'
                t = 0:0.0067*pi:2*pi;
          case '4P'
                t = 0:0.0067*pi:2*pi;
        end
    otherwise
                type = D;
                DS   = (type(1)*pi)/180;
                DE   = (type(2)*pi)/180;
                t    = DS:0.0067*pi:DE;
end
switch CF
    case 'clk'
                t = t;
    case 'aclk'
                t = sort(t,'descend');
end
for i = 1 : length(t)
r(i) = center(1) - radius*sin(t(i));
c(i) = center(2) - radius*cos(t(i));
% figure(1),
% plot(r,c,'b-');
end
varargout{1} = r;
varargout{2} = c;