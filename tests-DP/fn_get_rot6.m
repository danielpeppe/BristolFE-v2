function R6 = fn_get_rot6(th, about_axis)

% about_axis: x = 1, y = 2, z = 3

if about_axis == 1
    %TBC
elseif about_axis == 2
    R6 = [cosd(th)  0 sind(th) 0 0 0
          0         1 0        0 0 0
          -sind(th) 0 cosd(th) 0 0 0
          0         0          0 cosd(th)  0 sind(th)
          0         0          0 0         1 0
          0         0          0 -sind(th) 0 cosd(th)];
elseif about_axis == 3
    R6 = [cosd(th)  sind(th) 0 0 0 0
          -sind(th) cosd(th) 0 0 0 0
          0         0        1 0 0 0
          0         0        0 cosd(th)  sind(th) 0
          0         0        0 -sind(th) cosd(th) 0
          0         0        0 0         0        1];
else
    error('about_axis must equal 1, 2, or 3 (x, y, or z)')
end

end

