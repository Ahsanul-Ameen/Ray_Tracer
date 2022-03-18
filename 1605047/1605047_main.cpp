
#include "1605047_ccamera.h"
#include "1605047_classes.h"
#include "bitmap_image.h"


// declaration
int rec_level, dimensions;

vector <Object *> objects;
vector <Light *> lights;


void capture() {
	
	static int capture_no = 1;
	cout << "Capture [ " << capture_no << " ] " << endl;

	const int imageWidth = dimensions;
	const int imageHeight = dimensions;
	// initialize bitmap image and set background color
	bitmap_image image(imageWidth, imageHeight); 

	const double windowHeight = dimensions;
	const double windowWidth = dimensions;

	double planeDistance = (windowHeight/2.0) / tan(RAD(viewAngle/2.0));

	Point3D topLeft = eye + (l * planeDistance).endPoint; 
	topLeft	= topLeft - (r * (windowWidth / 2.0)).endPoint; 
	topLeft = topLeft + (u * (windowHeight / 2.0)).endPoint;


	double du = windowWidth / imageWidth;
	double dv = windowHeight / imageHeight;
 

	// Choose middle of the grid cell
	topLeft = topLeft + (r * (0.5 * du)).endPoint - (u * (0.5 * dv)).endPoint;


	int pixel_count = 0;

	int nearest;
	double t, t_min;


	for (int i = 0; i < imageWidth; ++i) {
		for (int j = 0; j < imageHeight; ++j) {
			image.set_pixel(i, j, 0, 0, 0);

			Point3D currentPixel = topLeft + (r * (i * du)).endPoint - (u * (j * dv)).endPoint;

			Vector3D direction(currentPixel - eye);
			direction.normalize();

			Ray ep_ray(eye, direction); 

			{
				double *color = new double[3];
				color[0] = color[1] = color[2] = 0.0;
				nearest = -1;
				t_min = INT_MAX;

				for(int k = 0, n = objects.size(); k < n; k++) {
					t = objects[k]->intersect(&ep_ray, color, 0); // level 0, dummy_color
					if(t > 0 && t < t_min) {
						t_min = t;
						nearest = k;
					}		
				}

				if(nearest != -1) {
					t_min = objects[nearest]->intersect(&ep_ray, color, 1); // level 1 ... rec_level

					pixel_count++;
					//cout << color[0] << " " << color[1] << " " << color[2] << endl;
					image.set_pixel(i, j, color[0] * 255.0, color[1] * 255.0, color[2] * 255.0);
				}


				delete[] color;
				color = nullptr;
			}
			
						
		}
	}	

    image.save_image("test.bmp");

    cout << "pixel_count = " << pixel_count << endl;
    cout << "image saved\n" << endl;
    ++capture_no;
}



int main(int argc, char **argv) {


	loadData(true);


    glutInit(&argc,argv);
    glutInitWindowSize(dimensions, dimensions);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("1605047_Offline_03");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occurring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}
