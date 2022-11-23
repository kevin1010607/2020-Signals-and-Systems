%% EECS2020 陳凱揚 108032053 Computer HW2 04/18/2021
% Function prototype of MyConv(), which implements convolution sum of x and h
%
%	Function input:
% 		x: input DT signal
%		support_x: support of x, which is a vector with 2 elements
%		h: impulse response of an LTI DT system
%		support_h: support of h, which a vector with 2 elements
%
%	Function output:
%		y: output DT signal of the LTI DT system
%		support_y: support of y, which is a vector with 2 elements	
%	
%	Note that you can add more input and output variables if needed
%	
%												Edited by Meng-Lin Li, 04/01/2021

function [y, support_y] = MyConv(x, support_x, h, support_h)
        % implement your own convolution sum here
        % the definition of "support", please see our lecture notes
        support_y = [0 0];
        support_y(1) = support_x(1)+support_h(1);
        support_y(2) = support_x(2)+support_h(2);
		x_len = length(x);
        h_len = length(h);
        y_len = x_len+h_len-1;
        y = zeros(1, y_len);
        for i = 1:x_len
            idx = i;
            for j = 1:h_len
                y(idx) = y(idx)+x(i)*h(j);
                idx = idx+1;
            end
        end
end	

