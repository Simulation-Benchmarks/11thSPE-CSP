Include "spe11a.geo";

lab_scale_bbox() = BoundingBox Surface{Surface{:}};
lab_scale_size_x = lab_scale_bbox(3) - lab_scale_bbox(0);
lab_scale_size_y = lab_scale_bbox(4) - lab_scale_bbox(1);

scale_factor_x = 8400.0/lab_scale_size_x;
scale_factor_y = 1200.0/lab_scale_size_y;

Dilate {{0.0, 0.0, 0.0}, {scale_factor_x, scale_factor_y, 1.0}}{Point{:};}
