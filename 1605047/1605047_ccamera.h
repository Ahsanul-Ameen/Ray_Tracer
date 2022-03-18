#ifdef _WIN32
    #include <windows.h>
#endif

#ifdef __APPLE__
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif


#include "1605047_classes.h"

#define PI (2.0 * acos(0.0))
#define RAD(d) (PI * d / 180.0)


#pragma once

extern void capture();
extern vector <Object *> objects;
extern vector <Light *> lights;



double viewAngle;
double mc, ac;
int draw_grid, draw_axes;
bool mouseClicked = false;


void drawAxes(int range = 100) {
    if(draw_axes == 1) {
        glColor3f(1.0, 0.2, 0.2);
        glBegin(GL_LINES); {
            glVertex3f( range,0,0);
            glVertex3f(-range,0,0);
        } glEnd();

        glColor3f(0.2, 1.0, 0.2);
        glBegin(GL_LINES); {
            glVertex3f(0,-range,0);
            glVertex3f(0, range,0);
        } glEnd();

        glColor3f(0.2, 0.2, 1.0);
        glBegin(GL_LINES); {
            glVertex3f(0,0, range);
            glVertex3f(0,0,-range);
        } glEnd();
    }
}



void drawGrid() {
    if(draw_grid == 1) {
        glColor3f(0.6, 0.6, 0.6);   //grey
        glBegin(GL_LINES); {
            for(int i =- 8; i <= 8; i++) {
                if(i==0)
                    continue;   //SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*10, -90, 0);
                glVertex3f(i*10,  90, 0);

                //lines parallel to X-axis
                glVertex3f(-90, i*10, 0);
                glVertex3f( 90, i*10, 0);
            }
        } glEnd();
    }
}

void drawSquare(double a) {
    glBegin(GL_QUADS); {
        glVertex3f(0, 0, -1);
        glVertex3f(0, +a, -1);
        glVertex3f(+a, +a, -1);
        glVertex3f(+a, 0, -1);
    } glEnd();
}

void drawCheckerBoard(const Point3D &bottom_left, double tileWidth) {
    if(draw_axes != 1) 
        return;
    int h = abs(2 * bottom_left.y / tileWidth) + 1e-7;
    int w = abs(2 * bottom_left.x / tileWidth) + 1e-7;

    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(((i + j) & 1) == 1) {
                glColor3f(1, 1, 1);
            } else {
                glColor3f(0, 0, 0);
            } 
            glPushMatrix(); {
                glTranslatef(bottom_left.x + j * tileWidth,
                         bottom_left.y + i * tileWidth, 
                         bottom_left.z);
                drawSquare(tileWidth);
            } glPopMatrix();
        }
    }
}

void drawTriangle(const vector<Point3D> &points) {
    glBegin(GL_TRIANGLES); {
        for(const Point3D & p : points) {
            glVertex3f(p.x, p.y, p.z);
        }
    } glEnd();    
}


void drawSphere(double radius, int slices, int stacks) {
    struct Point3D points[100][100];
    double h, r;
    //generate points
    for(int i = 0; i <= stacks; i++) {
        h = radius * sin(((double)i/(double)stacks) * (PI/2));
        r = radius * cos(((double)i/(double)stacks) * (PI/2));
        for(int j = 0; j <= slices; j++) {
            points[i][j].x = r * cos(((double)j/(double)slices) * 2*PI);
            points[i][j].y = r * sin(((double)j/(double)slices) * 2*PI);
            points[i][j].z = h;
        }
    }
    //draw quads using generated points
    for(int i = 0; i < stacks; i++) {
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
        for(int j = 0; j < slices; j++) {
            glBegin(GL_QUADS); {
                //upper hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
            } glEnd();
        }
    }
}


void keyboardListener(unsigned char key, int x, int y) {
    switch(key) {

        case '0':
            capture();
            break;
        
        case '1':   //  Rotate/Look left
            l = l * cos(RAD(ac)) + (u * l) * sin(RAD(ac));
            r = r * cos(RAD(ac)) + (u * r) * sin(RAD(ac));
            break;

        case '2':   // Rotate/Look right
            l = l * cos(-RAD(ac)) + (u * l) * sin(-RAD(ac));
            r = r * cos(-RAD(ac)) + (u * r) * sin(-RAD(ac));
            break;

        case '3':   // Look up
            u = u * cos(RAD(ac)) + (r * u) * sin(RAD(ac));
            l = l * cos(RAD(ac)) + (r * l) * sin(RAD(ac));
            break;

        case '4':   // Look down
            u = u * cos(-RAD(ac)) + (r * u) * sin(-RAD(ac));
            l = l * cos(-RAD(ac)) + (r * l) * sin(-RAD(ac));
            break;

        case '5':   // Tilt Clockwise
            u = u * cos(-RAD(ac)) + (l * u) * sin(-RAD(ac));
            r = r * cos(-RAD(ac)) + (l * r) * sin(-RAD(ac));
            break;

        case '6':   // Tilt Counterclockwise
            u = u * cos(RAD(ac)) + (l * u) * sin(RAD(ac));
            r = r * cos(RAD(ac)) + (l * r) * sin(RAD(ac));
            break;

        default:
            break;
    }
}


void specialKeyListener(int key, int x, int y) {
    switch(key) {
        case GLUT_KEY_UP:       // up arrow key
            eye = eye + (l*mc).endPoint;    // Move forward
            break;

        case GLUT_KEY_DOWN:     //down arrow key
            eye = eye - (l*mc).endPoint;    // Move backward
            break;

        case GLUT_KEY_LEFT:
            eye = eye - (r*mc).endPoint;    // Move left
            break;

        case GLUT_KEY_RIGHT:
            eye = eye + (r*mc).endPoint;    // Move right
            break;

        case GLUT_KEY_PAGE_UP:
            eye = eye + (u*mc).endPoint;    // Move up
            break;

        case GLUT_KEY_PAGE_DOWN:
            eye = eye - (u*mc).endPoint;    //  Move Down
            break;

        case GLUT_KEY_INSERT:
            break;

        case GLUT_KEY_HOME:
            break;
        
        case GLUT_KEY_END:
            break;

        default:
            break;
    }
}

void mouseListener(int button, int state, int x, int y) {   //x, y is the x-y of the screen (2D)
    switch(button) {
        case GLUT_LEFT_BUTTON:
            if(state == GLUT_DOWN) {        // 2 times?? in ONE click? -- solution is checking DOWN or UP
                mouseClicked = true;
            } else {
                mouseClicked = false;
            }
            break;

        case GLUT_RIGHT_BUTTON:
            if(state == GLUT_DOWN) {        // 2 times?? in ONE click? -- solution is checking DOWN or UP
                draw_axes = 1 - draw_axes;
            }
            break;

        case GLUT_MIDDLE_BUTTON:
            if(state == GLUT_DOWN) {        // 2 times?? in ONE click? -- solution is checking DOWN or UP
               // ...
            }
            break;

        default:
            break;
    }
}


void display() {

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);  //color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    Point3D lp = eye + l.endPoint;
    Point3D up = u.endPoint;
    gluLookAt(eye.x, eye.y, eye.z, lp.x, lp.y, lp.z, up.x, up.y, up.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/

    // draw objects
    for(const auto& o : objects) {
        o->draw();
    }

    // draw light sources
    for(const auto& l : lights) {
        l->draw();
    }

    // draw axes
    drawAxes(1000);

    drawGrid();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate() {
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init() {
    //codes for initialization
    draw_grid = 0;
    draw_axes = 1;

    eye.x = 100, eye.y = 100, eye.z = 0;
    u.endPoint.x = 0, u.endPoint.y = 0, u.endPoint.z = 1;
    r.endPoint.x = -1/sqrt(2.0), r.endPoint.y = 1/sqrt(2.0), r.endPoint.z = 0;
    l.endPoint.x = l.endPoint.y = -1/sqrt(2.0), l.endPoint.z = 0;

    viewAngle = 90.0; // fovY
    mc = 3, ac = 3; // degree

    mouseClicked = false;

    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(viewAngle,  1,  1,  1500.0);
    //field of view in the Y (vertically) in degree
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}