function R6_about_y = fn_get_rot6(th)

R6_about_y = [cosd(th)  0 sind(th) 0 0 0
              0         1 0        0 0 0
              -sind(th) 0 cosd(th) 0 0 0
              0         0          0 cosd(th)  0 sind(th)
              0         0          0 0         1 0
              0         0          0 -sind(th) 0 cosd(th)];

end

