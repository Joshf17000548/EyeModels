function  plotRay(ray)

rayStart = ray(:,1);
rayEnd =  ray(:,1) + 700.*ray(:,2);

ray = [rayStart, rayEnd]';

plot3(ray(:,1), ray(:,2), ray(:,3), '-');

end

