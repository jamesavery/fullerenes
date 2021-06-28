#pragma once;

typedef float real_t;
typedef uint16_t node_t; 


const int size = 60;
const int elements = 60*3;

//Data necessary to compute gradients for the C60ih Fullerene


const node_t cubic_neighbours_60[elements]  = {4,15,1,0,12,2,1,9,3,2,5,4,3,8,0,3,11,6,5,22,7,6,20,8,7,18,4,2,14,10,9,26,11,10,25,5,1,17,13,12,30,14,13,29,9,0,19,16,15,34,17,16,33,12,8,21,19,18,37,15,7,24,21,20,38,18,6,25,23,22,42,24,23,40,20,11,28,22,10,29,27,26,45,28,27,44,25,14,32,26,13,33,31,30,48,32,31,47,29,17,36,30,16,37,35,34,51,36,35,50,33,19,39,34,21,41,39,38,53,37,24,43,41,40,54,38,23,44,43,42,55,40,28,46,42,27,47,46,45,57,44,32,49,45,31,50,49,48,58,47,36,52,48,35,53,52,51,59,50,39,54,51,41,56,53,43,57,56,55,59,54,46,58,55,49,59,57,52,56,58};    
const node_t next_on_face_60[elements] = {8,16,2,15,13,3,12,10,4,9,6,0,5,18,1,2,25,7,11,23,8,22,21,4,20,19,3,1,29,11,14,27,5,26,22,3,0,33,14,17,31,9,30,26,2,4,37,17,19,35,12,34,30,1,7,38,15,21,34,0,6,40,18,24,39,8,5,28,24,25,43,20,42,41,7,10,44,6,9,32,28,29,46,25,45,42,11,13,47,10,12,36,32,33,49,29,48,45,14,16,50,13,15,39,36,37,52,33,51,48,17,18,53,16,20,54,37,41,51,19,23,55,38,43,53,21,22,46,40,44,56,24,27,57,23,26,49,44,47,55,28,31,58,27,30,52,47,50,57,32,35,59,31,34,54,50,53,58,36,38,56,35,40,59,39,42,58,54,57,52,41,45,59,43,48,56,46,51,55,49};
const node_t prev_on_face_60[elements] = {3,19,12,4,17,9,0,14,5,1,11,8,2,7,15,4,10,22,3,25,20,5,24,18,6,21,0,3,13,26,2,29,25,9,28,6,2,16,30,1,33,29,12,32,10,1,18,34,0,37,33,15,36,13,4,20,37,8,39,16,8,23,38,7,41,19,7,11,42,6,44,40,22,43,21,5,27,23,11,14,45,10,47,44,26,46,22,9,31,27,14,17,48,13,50,47,30,49,26,12,35,31,17,19,51,16,53,50,34,52,30,15,38,35,18,40,53,21,54,34,20,42,54,24,56,39,24,28,55,23,57,41,25,45,43,28,32,57,27,58,42,29,48,46,32,36,58,31,59,45,33,51,49,36,39,59,35,56,48,37,41,52,38,55,51,40,46,59,43,58,53,44,49,56,47,52,55,50,54,57};
const uint8_t face_right_60[elements]   = {6,6,5,6,6,5,6,6,5,6,6,5,6,6,5,6,5,6,5,6,6,6,5,6,5,6,6,6,5,6,5,6,6,6,5,6,6,5,6,5,6,6,6,5,6,6,5,6,5,6,6,6,5,6,5,6,6,6,5,6,6,6,5,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,6,5,6,5,6,6,5,6,6,6,5,6,6,6,5,6,6,5,5,6,6,6,6,5,6,6,5,5,6,6,6,6,5,6,6,5,5,6,6,6,6,5,6,6,5,5,6,6,6,6,5,6,5,6,5,6,6,6,5,6,6,5,6,6,5,6};
const real_t X_60[elements] = {2.82489,8.43261e-17,14.2017,8.0219,0.604118,12.0396,10.4266,6.01981,8.04461,4.53413,6.64511,12.0396,1.41245,2.44643,14.2017,3.29447,11.5801,8.04461,-3.29447,11.5801,8.04461,-4.53413,6.64511,12.0396,-1.41245,2.44643,14.2017,12.299,7.10085,2.82489,9.58483,10.4795,-2.82489,4.63059,13.4256,2.82489,11.6759,-2.93695,8.04461,13.9422,-2.70257,2.82489,13.8679,3.06098,-2.82489,1.41245,-2.44643,14.2017,3.48777,-7.24923,12.0396,8.38143,-8.64315,8.04461,-2.82489,5.87089e-16,14.2017,-1.41245,-2.44643,14.2017,-10.4266,6.01981,8.04461,-8.0219,0.604118,12.0396,-4.63059,13.4256,2.82489,-9.58483,10.4795,-2.82489,-12.299,7.10085,2.82489,2.39962e-16,14.2017,-2.82489,8.97979,8.0197,-8.04461,3.58009,7.20408,-12.0396,1.07573e-16,12.0396,-8.04461,11.4352,3.76688,-8.04461,12.299,-7.10085,-2.82489,10.4266,-6.01981,-8.04461,8.02896,-0.501593,-12.0396,9.31158,-10.723,2.82489,-2.56576e-15,-12.0396,8.04461,-2.91345e-15,-14.2017,2.82489,4.28306,-13.5404,-2.82489,-3.48777,-7.24923,12.0396,-11.6759,-2.93695,8.04461,-8.38143,-8.64315,8.04461,-13.8679,3.06098,-2.82489,-13.9422,-2.70257,2.82489,-8.97979,8.0197,-8.04461,-11.4352,3.76688,-8.04461,-3.58009,7.20408,-12.0396,1.33536,2.48934,-14.2017,-1.33536,2.48934,-14.2017,2.82351,-0.088213,-14.2017,4.44887,-6.70249,-12.0396,1.48815,-2.40113,-14.2017,2.45537,-11.7866,-8.04461,-4.28306,-13.5404,-2.82489,-2.45537,-11.7866,-8.04461,-9.31158,-10.723,2.82489,-12.299,-7.10085,-2.82489,-8.02896,-0.501593,-12.0396,-10.4266,-6.01981,-8.04461,-2.82351,-0.088213,-14.2017,-1.48815,-2.40113,-14.2017,-4.44887,-6.70249,-12.0396};
