// Skeleton code
#include <windows.h>
#include <time.h>
#include <math.h>
#define GLUT_DISABLE_ATEXIT_HACK // In case of errors like undefined reference to `__glutInitWithExit@12’, should be before including glut.h 
#include <gl\gl.h>
#include <gl\glu.h>
#include <gl\glut.h>
#include <fstream>
#include <iostream>

#define PI 3.14159265
int screenWidth = 600;
int screenHeight = 600;
int delay = 10;

double alfa = 0;
int    fx = 0, fy = 0, fz = 0;
float  sphi = 0.0, stheta = 0.0;
float  sside = 0, sdepth = -5;
float  sx = 0, sy = 0;
bool  mouse_left_click, mouse_middle_click, mouse_right_click;
int   mouseX, mouseY;

// -------------------------- radius ------------------------------------
double r = 1;

// -------------------------- theta and tangle  --------------------------
double theta = 0;
double tangle = 0;

// --------------------- LOOK AT SETTINGS -------------------
double lookat[9] = { 0,-25,-18, 0, 0, 0, 0, 1, 0 };

// ------------ contastant k --------------
// You can change k to change the entire shape size. When the window is opened, you can enlarge or reduce it by right clicking.
double k = 1;


// ----------------------------- Tires ------------------------------------
// ---------------- tire_front-left ----------------
double t_fl_1[3] = { 0,0,0 };
double t_fl_2[3] = { 0,0,0 };
// The vertex we will use to make the line turn left and right.
double t_fl_11[3] = { 0,0,0 };
double t_fl_22[3] = { 0,0,0 };
// ---------------- tire_front-right ----------------
double t_fr_1[3] = { 0,0,0 };
double t_fr_2[3] = { 0,0,0 };
// The vertex we will use to make the line turn left and right.
double t_fr_11[3] = { 0,0,0 };
double t_fr_22[3] = { 0,0,0 };
// ---------------- tire_back-left ----------------
double t_bl_1[3] = { 0,0,0 };
double t_bl_2[3] = { 0,0,0 };
// The vertex we will use to make the line turn left and right.
double t_bl_11[3] = { 0,0,0 };
double t_bl_22[3] = { 0,0,0 };
// ---------------- tire_back-right ----------------
double t_br_1[3] = { 0,0,0 };
double t_br_2[3] = { 0,0,0 };
// The vertex we will use to make the line turn left and right.
double t_br_11[3] = { 0,0,0 };
double t_br_22[3] = { 0,0,0 };



// _______________________________________________________________________________________            _______________________________________________________________________________________ 
// _______________________________________________________________________________________  STARTÝNG  _______________________________________________________________________________________          


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<< tires_animation >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//   Animation routine which calls itself after “delay” miliseconds.
void tires_animation(int frame)
{

    alfa += 10;
    if (alfa > 360) alfa -= 360;

    // ------------------------------------------  tire_front-left  ----------------------------------------------------------------------
    // The necessary formula to give the impression of moving.
    t_fl_1[0] = 5 * k ;
    t_fl_1[1] = -1 * k + k * cos(alfa * PI / 180);
    t_fl_1[2] = 2.1 * k + k * sin(alfa * PI / 180);

    t_fl_2[0] = 5 * k;
    t_fl_2[1] = -1 * k + k * cos((alfa + 180) * PI / 180);
    t_fl_2[2] = 2.1 * k + k * sin((alfa + 180) * PI / 180);

    // The back tires is depend to the 'theta' as it turn left and right. Rotate it 45 degrees to the left, as its normal state is not what I want.
    t_fl_11[0] = 5 * k + (t_fl_1[0] - 5 * k) * cos((theta - 45) * PI / 180) - (t_fl_1[1] - (-1 * k)) * sin((theta - 45) * PI / 180);
    t_fl_11[1] = (-1 * k) + (t_fl_1[0] - 5 * k) * sin((theta - 45) * PI / 180) + (t_fl_1[1] - (-1 * k)) * cos((theta - 45) * PI / 180);
    t_fl_11[2] = t_fl_1[2];

    t_fl_22[0] = 5 * k + (t_fl_2[0] - 5 * k) * cos((theta - 45) * PI / 180) - (t_fl_2[1] - (-1 * k)) * sin((theta - 45) * PI / 180);
    t_fl_22[1] = (-1 * k) + (t_fl_2[0] - 5 * k) * sin((theta - 45) * PI / 180) + (t_fl_2[1] - (-1 * k)) * cos((theta - 45) * PI / 180);
    t_fl_22[2] = t_fl_2[2];
    // ------------------------------------------------------------------------------------------------------------------------------------


    // ------------------------------------------  tire_front-right  ----------------------------------------------------------------------
    // The necessary formula to give the impression of moving.
    t_fr_1[0] = -1 * k ;
    t_fr_1[1] = 5 * k + k * cos(alfa * PI / 180);
    t_fr_1[2] = 2.1 * k + k * sin(alfa * PI / 180);

    t_fr_2[0] = -1 * k ;
    t_fr_2[1] = 5 * k + k * cos((alfa + 180) * PI / 180);
    t_fr_2[2] = 2.1 * k + k * sin((alfa + 180) * PI / 180);

    // The back tires is depend to the 'theta' as it turn left and right. Rotate it 45 degrees to the left, as its normal state is not what I want.
    t_fr_11[0] = (-1 * k) + (t_fr_1[0] - (-1 * k)) * cos((theta - 45) * PI / 180) - (t_fr_1[1] - (5 * k)) * sin((theta - 45) * PI / 180);
    t_fr_11[1] = (5 * k) + (t_fr_1[0] - (-1 * k)) * sin((theta - 45) * PI / 180) + (t_fr_1[1] - (5 * k)) * cos((theta - 45) * PI / 180);
    t_fr_11[2] = t_fr_1[2];

    t_fr_22[0] = (-1 * k) + (t_fr_2[0] - (-1 * k)) * cos((theta - 45) * PI / 180) - (t_fr_2[1] - (5 * k)) * sin((theta - 45) * PI / 180);
    t_fr_22[1] = (5 * k) + (t_fr_2[0] - (-1 * k)) * sin((theta - 45) * PI / 180) + (t_fr_2[1] - (5 * k)) * cos((theta - 45) * PI / 180);
    t_fr_22[2] = t_fr_2[2];
    // -------------------------------------------------------------------------------------------------------------------------------------


    // ------------------------------------------  tire_back-left  -------------------------------------------------------------------------
    // The necessary formula to give the impression of moving.
    t_bl_1[0] = 1 * k ;
    t_bl_1[1] = -5 * k + k * cos(alfa * PI / 180);
    t_bl_1[2] = 2.1 * k + k * sin(alfa*PI / 180);

    t_bl_2[0] = 1 * k ;
    t_bl_2[1] = -5 * k + k * cos((alfa + 180) * PI / 180);
    t_bl_2[2] = 2.1 * k + k * sin((alfa + 180) * PI / 180);
    // The back tires is not depend to the 'theta' as it does not turn left and right. Rotate it 45 degrees to the left, as its normal state is not what I want.
    t_bl_11[0] = (1 * k) + (t_bl_1[0] - (1 * k)) * cos((- 45) * PI / 180) - (t_bl_1[1] - (-5 * k)) * sin((- 45) * PI / 180);
    t_bl_11[1] = (-5 * k) + (t_bl_1[0] - (1 * k)) * sin((- 45) * PI / 180) + (t_bl_1[1] - (-5 * k)) * cos((- 45) * PI / 180);
    t_bl_11[2] = t_bl_1[2];

    t_bl_22[0] = (1 * k) + (t_bl_2[0] - (1 * k)) * cos((- 45) * PI / 180) - (t_bl_2[1] - (-5 * k)) * sin((- 45) * PI / 180);
    t_bl_22[1] = (-5 * k) + (t_bl_2[0] - (1 * k)) * sin((- 45) * PI / 180) + (t_bl_2[1] - (-5 * k)) * cos((- 45) * PI / 180);
    t_bl_22[2] = t_bl_2[2];
    // --------------------------------------------------------------------------------------------------------------------------------------


    // -------------------------------------------  tire_back-right  ------------------------------------------------------------------------
    // The necessary formula to give the impression of moving.
    t_br_1[0] = -5 * k;
    t_br_1[1] = k + k * cos(alfa * PI / 180);
    t_br_1[2] = 2.1 * k + k * sin(alfa * PI / 180);

    t_br_2[0] = -5 * k;
    t_br_2[1] = k + k * cos((alfa + 180) * PI / 180);
    t_br_2[2] = 2.1 * k + k * sin((alfa + 180) * PI / 180);
    // The back tires is not depend to the 'theta' as it does not turn left and right. Rotate it 45 degrees to the left, as its normal state is not what I want.
    t_br_11[0] = (-5 * k) + (t_br_1[0] - (-5 * k)) * cos((-45) * PI / 180) - (t_br_1[1] - (k)) * sin((-45) * PI / 180);
    t_br_11[1] = (k) + (t_br_1[0] - (-5 * k)) * sin((-45) * PI / 180) + (t_br_1[1] - (k)) * cos((-45) * PI / 180);
    t_br_11[2] = t_br_1[2];

    t_br_22[0] = (-5 * k) + (t_br_2[0] - (-5 * k)) * cos((-45) * PI / 180) - (t_br_2[1] - (k)) * sin((-45) * PI / 180);
    t_br_22[1] = (k) + (t_br_2[0] - (-5 * k)) * sin((-45) * PI / 180) + (t_br_2[1] - (k)) * cos((-45) * PI / 180);
    t_br_22[2] = t_br_2[2];
    // --------------------------------------------------------------------------------------------------------------------------------------

    // Calling Itself
    glutTimerFunc(delay, tires_animation, 0);
    glutPostRedisplay();
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  circle drawing function for front tires >>>>>>>>>>>>>>>>>>>>>>>>>>>
void DrawCircleForFrontTires(double* center)
{
    int num_segments = 20;
    // center points
    double cx = center[0];
    double cy = center[1];
    double cz = center[2];
    float tbeta = 2 * PI / float(num_segments);
    float c = cosf(tbeta);//precalculate the sine and cosine
    float s = sinf(tbeta);
    float t;

    float x = r*k;//we start at angle = 0 
    float y = 0;
    double array_x[21];
    double array_y[21];
    double array_z[21];
    // We take the points around the circle and throw them into arrays.
    for (int i = 0; i < num_segments; i++)
    {
        array_x[i] = x + cx;
        array_y[i] = x + cy;
        array_z[i] = y + cz;

        //apply the rotation matrix
        t = x;
        x = x * c - y * s;
        y = t * s + y * c;
    }
    // draw the circle using gl_lines
    glBegin(GL_LINES);

    double vektor1[3] = { 0,0,0 };
    double vektor2[3] = { 0,0,0 };
    
    // draw the circle using arrays
    for (int j = 0; j < 20; j++)
    { 
        // We use the following formula
        // x' = x.cosQ - y.sinQ
        // y' = x.sinQ + y.cosQ
        // z' = z
        if (j == 19) { // We draw a line between the last point and the first point to draw the last line not drawn.
            
            vektor1[0] = cx + (array_x[19] - cx) * cos(theta * PI / 180) - (array_y[19] - cy) * sin(theta * PI / 180);
            vektor1[1] = cy + (array_x[19] - cx) * sin(theta * PI / 180) + (array_y[19] - cy) * cos(theta * PI / 180);
            vektor1[2] = array_z[19];

            vektor2[0] = cx + (array_x[0] - cx) * cos(theta * PI / 180) - (array_y[0] - cy) * sin(theta * PI / 180);
            vektor2[1] = cy + (array_x[0] - cx) * sin(theta * PI / 180) + (array_y[0] - cy) * cos(theta * PI / 180);
            vektor2[2] = array_z[0];


            glVertex3dv(vektor1);
            glVertex3dv(vektor2);
        }
        else {

            vektor1[0] = cx + (array_x[j] - cx) * cos(theta * PI / 180) - (array_y[j] - cy) * sin(theta * PI / 180);
            vektor1[1] = cy + (array_x[j] - cx) * sin(theta * PI / 180) + (array_y[j] - cy) * cos(theta * PI / 180);
            vektor1[2] = array_z[j];

            vektor2[0] = cx + (array_x[j + 1] - cx) * cos(theta * PI / 180) - (array_y[j + 1] - cy) * sin(theta * PI / 180);
            vektor2[1] = cy + (array_x[j + 1] - cx) * sin(theta * PI / 180) + (array_y[j + 1] - cy) * cos(theta * PI / 180);
            vektor2[2] = array_z[j + 1];

            glVertex3dv(vektor1);
            glVertex3dv(vektor2);
        }
    }
    glEnd();

}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< circle drawing function for back tires >>>>>>>>>>>>>>>>>>>>>>>>>>>
void DrawCircleForBackTires(double* center)
{
    int num_segments = 20;
    // center points
    double cx = center[0];
    double cy = center[1];
    double cz = center[2];
    float tbeta = 2 * PI / float(num_segments);
    float c = cosf(tbeta);//precalculate the sine and cosine
    float s = sinf(tbeta);
    float t;

    float x = r*k;//we start at angle = 0 
    float y = 0;
    double array_x[21];
    double array_y[21];
    double array_z[21];

    // We take the points around the circle and throw them into arrays.
    for (int i = 0; i < num_segments; i++)
    {
        array_x[i] = x + cx;
        array_y[i] = x + cy;
        array_z[i] = y + cz;

        //apply the rotation matrix
        t = x;
        x = x * c - y * s;
        y = t * s + y * c;
    }

    // draw the circle using gl_lines
    glBegin(GL_LINES);

    double vektor1[3] = { 0,0,0 };
    double vektor2[3] = { 0,0,0 };

    // draw the circle using arrays
    for (int j = 0; j < 20; j++)
    {
        if (j == 19) {
            vektor1[0] = array_x[19];
            vektor1[1] = array_y[19];
            vektor1[2] = array_z[19];

            vektor2[0] = array_x[0];
            vektor2[1] = array_y[0];
            vektor2[2] = array_z[0];

            glVertex3dv(vektor1);
            glVertex3dv(vektor2);
        }
        else {
            vektor1[0] = array_x[j];
            vektor1[1] = array_y[j];
            vektor1[2] = array_z[j];

            vektor2[0] = array_x[j + 1];
            vektor2[1] = array_y[j + 1];
            vektor2[2] = array_z[j + 1];

            glVertex3dv(vektor1);
            glVertex3dv(vektor2);
        }
    }
    glEnd();
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< myinit >>>>>>>>>>>>>>>>>>>>>>>>>>>
void myInit(double* elavation)
{
    glColor3f(0.0f, 0.0f, 0.0f);  // set color of stuff
    glShadeModel(GL_FLAT);	// or can be GL_SMOOTH

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // Produce the perspective projection
    gluPerspective(45.0f, 1.0, 1.0, 100.0);
    // Car view menu needs these arrays. Because we will change If we want to change the view of the car.
    gluLookAt(elavation[0], elavation[1], elavation[2], elavation[3], elavation[4], elavation[5], elavation[6], elavation[7], elavation[8]);

    glMatrixMode(GL_MODELVIEW);

    // Start animation
    tires_animation(0);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Tires lines vertices for drawing >>>>>>>>>>>>>>>>>>>>>>>>>>>
void tires_lines() {
    // Tires lines vertices for drawing
    glVertex3dv(t_fl_11);
    glVertex3dv(t_fl_22);

    glVertex3dv(t_fr_11);
    glVertex3dv(t_fr_22);

    glVertex3dv(t_bl_11);
    glVertex3dv(t_bl_22);

    glVertex3dv(t_br_11);
    glVertex3dv(t_br_22);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< tires circle drawing and rotation >>>>>>>>>>>>>>>>>>>>>>>>>>>
void tires_circle() {
    // Tires circle center points
    double tire1[3] = { 5 * k, -1 * k, 2.1 * k};
    double tire2[3] = { -1 * k, 5* k, 2.1 * k };
    double tire3[3] = { k, -5 * k, 2.1 * k };
    double tire4[3] = { -5 * k, k, 2.1 * k };

    // The front tires will turn
    DrawCircleForFrontTires(tire1);
    DrawCircleForFrontTires(tire2);
    // The back tires will not turn
    DrawCircleForBackTires(tire3);
    DrawCircleForBackTires(tire4);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< steering wheel vertices for rotation and drawing >>>>>>>>>>>>>>>>>>>>>>>>>>>
void steering_wheel() {
    // -------------------------- Steering Wheel ------------------------------
// SteeringWheel_left
    double sw_l[3] = { 3 * k, 2 * k, -3.1 * k };
    // SteeringWheel_right
    double sw_r[3] = { 1 * k, 2 * k, -3.1 * k };
    // SteeringWheel_bottom
    double sw_b[3] = { 2 * k, 1 * k, -3.1 * k };
    // SteeringWheel_top
    double sw_t[3] = { 2 * k, 3 * k, -3.1 * k };

    // Steering wheel rotation (depend on theta)
    // Changed the steering lines from 'x' to '+' shape by changing angles. We rotate the starting position 45 degrees to the left.
    sw_l[0] = 2 * k + k * (cos((theta - 45) * PI / 180));
    sw_l[1] = 2 * k + k * (sin((theta - 45) * PI / 180));

    sw_r[0] = 2 * k + k * (cos((theta + 135) * PI / 180));
    sw_r[1] = 2 * k + k * (sin((theta + 135) * PI / 180));

    sw_b[0] = 2 * k + k * (cos((theta - 135) * PI / 180));
    sw_b[1] = 2 * k + k * (sin((theta - 135) * PI / 180));

    sw_t[0] = 2 * k + k * (cos((theta + 45) * PI / 180));
    sw_t[1] = 2 * k + k * (sin((theta + 45) * PI / 180));

    // steering wheel lines vertices for drawing
    glVertex3dv(sw_l);
    glVertex3dv(sw_r);

    glVertex3dv(sw_b);
    glVertex3dv(sw_t);

    // first point of circle for drawing
    // z point will not change (z = -3.1 * k)
    double x = 1;
    double y = 2;
    double inc = 0.05; //  Amount of increase and decrease

    // <-----------------  My own algorithm to draw a circle using formula  ----------------------->
    //                     -----------------------------------------------
    // Find vertices for drawing Circle using this formula  = (x-2)^2 + (y-2)^2 = 1   ::> Cricle Equation     cx = 2 , cy = 2 , r = 1

    // half of the circle(lower part)         formula for lower part = 2 - sqrt(-3 + 4 x - x^2)
    for (int i = 1; i <= 40; i++) {
        // Find vertices for drawing the first quarter of the circle
        if (i <= 20) {
            double c1[3] = { x * k,y * k,-3.1 * k };
            x = x + inc;
            if (x == 2) {
                y = 1;
            }
            else {
                y = (2 - sqrt(-pow(x, 2) + 4 * x - 3));
            }
            double c2[3] = { x * k,y * k,-3.1 * k };

            glVertex3dv(c1);
            glVertex3dv(c2);
        }
        // Find vertices for drawing the second quarter of the circle
        else {
            double c1[3] = { x * k,y * k,-3.1 * k };
            x = x + inc;
            if (x == 3) {
                y = 2;
            }
            else {
                y = (2 - sqrt(-pow(x, 2) + 4 * x - 3));
            }
            double c2[3] = { x * k,y * k,-3.1 * k };

            glVertex3dv(c1);
            glVertex3dv(c2);
        }
    }
    // half of the circle (upper part)      formula for upper part = 2 + sqrt(-3 + 4 x - x^2)
    for (int j = 1; j <= 40; j++) {
        // Find vertices for drawing the third quarter of the circle
        if (j <= 20) {
            double c1[3] = { x * k,y * k,-3.1 * k };
            x = x - inc;
            if (x == 2) {
                y = 3;
            }
            else {
                y = (2 + sqrt(-pow(x, 2) + 4 * x - 3));
            }
            double c2[3] = { x * k,y * k,-3.1 * k };

            glVertex3dv(c1);
            glVertex3dv(c2);
        }
        // Find vertices for drawing the fourth quarter of the circle
        else {
            double c1[3] = { x * k,y * k,-3.1 * k };
            x = x - inc;
            if (x == 1) {
                y = 2;
            }
            else {
                y = (2 + sqrt(-pow(x, 2) + 4 * x - 3));
            }
            double c2[3] = { x * k,y * k,-3.1 * k };

            glVertex3dv(c1);
            glVertex3dv(c2);
        }
    }

}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< steering column vertices for drawing >>>>>>>>>>>>>>>>>>>>>>>>>>>
void steering_column() {

    // -------------------------- Steering Column -----------------------------
    // SteeringColumn_bottom
    double sc_b[3] = { 2 * k, 2 * k,  2.1 * k };
    // SteeringColumn_top
    double sc_t[3] = { 2 * k, 2 * k, -3.1 * k };

    // Steering column vertices for drawing
    glVertex3dv(sc_b);
    glVertex3dv(sc_t);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< axles vertices for drawing >>>>>>>>>>>>>>>>>>>>>>>>>>>
void axles() {
    // ------------------------------- Axles ----------------------------------
// axles_front-left
    double axles_fl[3] = { 5 * k, -k, 2.1 * k };
    // axles_front-right
    double axles_fr[3] = { -k, 5 * k, 2.1 * k };
    // axles_back-left
    double axles_bl[3] = { k, -5 * k, 2.1 * k };
    // axles_back-right
    double axles_br[3] = { -5 * k, k, 2.1 * k };

    // front axles rotation
    axles_fl[0] = 2 * k + 3 * sqrt(2) * (k * (cos((tangle - 45) * PI / 180)));
    axles_fl[1] = 2 * k + 3 * sqrt(2) * (k * (sin((tangle - 45) * PI / 180)));

    axles_fr[0] = 2 * k + 3 * sqrt(2) * (k * (cos((tangle + 135) * PI / 180)));
    axles_fr[1] = 2 * k + 3 * sqrt(2) * (k * (sin((tangle + 135) * PI / 180)));



    // ------- axle between front wheels --- vertices
    glVertex3dv(axles_fl);
    glVertex3dv(axles_fr);
    // ------- axle between back wheels --- vertices
    glVertex3dv(axles_bl);
    glVertex3dv(axles_br);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< cube vertices for drawing >>>>>>>>>>>>>>>>>>>>>>>>>>>
void cube() {
    // --------------------------------- Cube ------------------------------------
// --- Front surface of the cube ---

// cube_front-left-top
    double cube_flt[3] = { 4 * k, 0, 0 };
    // cube_front-right-top
    double cube_frt[3] = { 0, 4 * k, 0 };
    // cube_front-right-bottom
    double cube_frb[3] = { -4 * k, 0, 0 };
    // cube_front-left-bottom
    double cube_flb[3] = { 0, -4 * k, 0 };

    // --- The back surface of the cube ---

    // cube_back-left-top
    double cube_blt[3] = { 4 * k, 0, 2 * k };
    // cube_back-right-top
    double cube_brt[3] = { 0, 4 * k, 2 * k };
    // cube_back-right-bottom
    double cube_brb[3] = { -4 * k, 0, 2 * k };
    // cube_back-left-bottom
    double cube_blb[3] = { 0, -4 * k, 2 * k };
    // ------- Cube front face -------------- vertices
    glVertex3dv(cube_flt);
    glVertex3dv(cube_frt);
    glVertex3dv(cube_frt);
    glVertex3dv(cube_frb);
    glVertex3dv(cube_frb);
    glVertex3dv(cube_flb);
    glVertex3dv(cube_flb);
    glVertex3dv(cube_flt);
    // --------- Cube back face --------------- vertices
    glVertex3dv(cube_blt);
    glVertex3dv(cube_brt);
    glVertex3dv(cube_brt);
    glVertex3dv(cube_brb);
    glVertex3dv(cube_brb);
    glVertex3dv(cube_blb);
    glVertex3dv(cube_blb);
    glVertex3dv(cube_blt);
    // --------- Cube connecting lines between the front surface and the back surface ----- vertices
    glVertex3dv(cube_flt);
    glVertex3dv(cube_blt);
    glVertex3dv(cube_frt);
    glVertex3dv(cube_brt);
    glVertex3dv(cube_frb);
    glVertex3dv(cube_brb);
    glVertex3dv(cube_flb);
    glVertex3dv(cube_blb);
}

//<<<<<<<<<<<<<<<<<<<<<<<< myKeyboard >>>>>>>>>>>>>>>>>>>>>>
void myKeyboard(unsigned char key, int x, int y)
{
    switch (key) {

    case 27:  // Escape
        exit(-1);
    case 'a': // If 'a' is pressed it will turn right
        theta -= 1; // Required angle to turn left and right
        tangle -= 0.12; // Required angle to axles
        if (theta < -360) theta += 360;
        if (tangle < -360) tangle += 360;
        break;
    case 's': // If 'a' is pressed it will turn left
        theta += 1;
        tangle += 0.12;
        if (theta > 360) theta -= 360;
        if (tangle > 360) tangle -= 360;
        break;

    case 'A': // If 'A' is pressed it will turn right
        theta -= 1;
        tangle -= 0.12;
        if (theta < -360) theta += 360;
        if (tangle < -360) tangle += 360;
        break;
    case 'S': // If 'S' is pressed it will turn left
        theta += 1;
        tangle += 0.12;
        if (theta > 360) theta -= 360;
        if (tangle > 360) tangle -= 360;
        break;
    }
  
    glutPostRedisplay();
}

//<<<<<<<<<<<<<<<<<<<<<<<<<< myReshape >>>>>>>>>>>>>>>>>>>
void myReshape(int width, int height)
{ // adjust the camera aspect ratio to match that of the viewport
    glViewport(0, 0, width, height); // update viewport
    //glOrtho(-width,width,-height,height,-1000,1000);
    glOrtho(-1, 1, -1, 1, -1, 1);

}

//<<<<<<<<<<<<<<<<<<<<<<<<<< drawText >>>>>>>>>>>>>>>>>>>
void drawText(const char* text, int lenght, int x, int y, int z) {

    glMatrixMode(GL_PROJECTION);
    double* matrix = new double[200];
    glGetDoublev(GL_PROJECTION_MATRIX, matrix);
    glLoadIdentity();
    glOrtho(0, 800, 0, 600, -5, 5);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glLoadIdentity();
    glRasterPos3i(x, y, z);
    for (int i = 0; i < lenght; i++) {

        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, (int)text[i]);
    }
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd(matrix);
    glMatrixMode(GL_MODELVIEW);

}

//<<<<<<<<<<<<<<<<<<<<<<< myDisplay >>>>>>>>>>>>>>>>>>>>>>>>>>
void myDisplay(void)
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// clear screen

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glTranslatef(sside, 0, -sdepth);
    glRotatef(-stheta, 1, 0, 0);
    glRotatef(sphi, 0, 1, 0);
    glTranslatef(sx, 0, -sy);

    glColor3f(1, 1, 1);

    glBegin(GL_LINES);

    // my functions for drawing
    tires_lines();
    cube();
    axles();
    steering_column();
    steering_wheel();

    glEnd();

    tires_circle();

    // draw text
    std::string text1;
    text1 = "Right click to change the view and size";
    drawText(text1.data(), text1.size(), 150, 450, 0);

 
    glutSwapBuffers();
}

//<<<<<<<<<<<<<<<<<<<<<<< myDisplayForSecondWindow >>>>>>>>>>>>>>>>>>>>>>>>>>
void myDisplayForSecondWindow() {
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// clear screen

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glColor3f(1, 1, 1);

    glBegin(GL_LINES);
    // my functions for drawing
    tires_lines();
    cube();
    axles();
    steering_column();
    steering_wheel();
    glEnd();

    tires_circle();

    glutSwapBuffers();
}

/** myMouse()
 *
 * This event callback is executed whenever there is a mouse event
 */
void myMouse(int button, int state, int x, int y)
{
    mouseX = x; mouseY = y;
    mouse_left_click = ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN));
    mouse_middle_click = ((button == GLUT_MIDDLE_BUTTON) &&
        (state == GLUT_DOWN));
    mouse_right_click = ((button == GLUT_RIGHT_BUTTON) &&
        (state == GLUT_DOWN));
    glutPostRedisplay();
}


/** myMouseMove
 *
 * This even callback is executed whenver the mouse is moved
 */
void myMouseMove(int x, int y) {
    // rotate
    if (mouse_left_click)
    {
        sphi += (float)(x - mouseX) / 4.0;
        stheta += (float)(mouseY - y) / 4.0;
        // if (stheta<0) stheta=0;
    }
    mouseX = x;
    mouseY = y;
    glutPostRedisplay();
}

//<<<<<<<<<<<<<<<<<<<<<<< Car View Menu  >>>>>>>>>>>>>>>>>>>>>>>>>>
void CarViewMenu(int value) {

    switch (value)
    {
        // case 1 : We change the "gluLookat" for the top view of the car.   
    case 1:
        lookat[0] = 0;
        lookat[1] = 0;
        lookat[2] = -25;
        lookat[3] = 0;
        lookat[4] = 0;
        lookat[5] = 0;
        lookat[6] = 0;
        lookat[7] = 1;
        lookat[8] = 0;

        // New window and window position
        glutInitWindowPosition(700, 50);
        glutCreateWindow("TOP");

        glutKeyboardFunc(myKeyboard);
        glutReshapeFunc(myReshape);
        glutDisplayFunc(myDisplayForSecondWindow);
        myInit(lookat);

        break;
        // case 2 : We change the "gluLookat" for the bottom view of the car.
    case 2:
        lookat[0] = 0;
        lookat[1] = 0;
        lookat[2] = 30;
        lookat[3] = 0;
        lookat[4] = 0;
        lookat[5] = 0;
        lookat[6] = 0;
        lookat[7] = 1;
        lookat[8] = 0;

        // New window and window position
        glutInitWindowPosition(700, 50);
        glutCreateWindow("BOTTOM");

        glutKeyboardFunc(myKeyboard);
        glutReshapeFunc(myReshape);
        glutDisplayFunc(myDisplayForSecondWindow);
        myInit(lookat);
        break;
        // case 3 : We change the "gluLookat" for the right view of the car.
    case 3:
        lookat[0] = -25;
        lookat[1] = 25;
        lookat[2] = 0;
        lookat[3] = 0;
        lookat[4] = 0;
        lookat[5] = 0;
        lookat[6] = 0;
        lookat[7] = 0;
        lookat[8] = -1;

        // New window and window position
        glutInitWindowPosition(700, 50);
        glutCreateWindow("RIGHT");

        glutKeyboardFunc(myKeyboard);
        glutReshapeFunc(myReshape);
        glutDisplayFunc(myDisplayForSecondWindow);
        myInit(lookat);
        break;
        // case 4 : We change the "gluLookat" for the left view of the car.
    case 4:
        lookat[0] = 25;
        lookat[1] = -25;
        lookat[2] = 0;
        lookat[3] = 0;
        lookat[4] = 0;
        lookat[5] = 0;
        lookat[6] = 0;
        lookat[7] = 0;
        lookat[8] = -1;

        // New window and window position
        glutInitWindowPosition(700, 50);
        glutCreateWindow("LEFT");

        glutKeyboardFunc(myKeyboard);
        glutReshapeFunc(myReshape);
        glutDisplayFunc(myDisplayForSecondWindow);
        myInit(lookat);

        break;
        // case 5 : We change the "gluLookat" for the normal view of the car.
    case 5:
        lookat[0] = 0;
        lookat[1] = -25;
        lookat[2] = -18;
        lookat[3] = 0;
        lookat[4] = 0;
        lookat[5] = 0;
        lookat[6] = 0;
        lookat[7] = 0;
        lookat[8] = -1;

        k = 1;
        // New window and window position
        glutInitWindowPosition(700, 50);
        glutCreateWindow("NORMALLY");

        glutKeyboardFunc(myKeyboard);
        glutReshapeFunc(myReshape);
        glutDisplayFunc(myDisplayForSecondWindow);
        myInit(lookat);
        break;
        // case 6 : We change the "k" for the larger view of the car.
    case 6:
        k = k + 0.2;
        break;
        // case 7: We change the "k" for the smaller view of the car.
    case 7:
        k = k - 0.2;
        break;
        
    }
    glutPostRedisplay();
}

//<<<<<<<<<<<<<<<<<<<<<<< main >>>>>>>>>>>>>>>>>>>>>>>>>>
void main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(screenWidth, screenHeight);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("WireCar");
    glutKeyboardFunc(myKeyboard);
    glutReshapeFunc(myReshape);
    glutDisplayFunc(myDisplay);
    glutMouseFunc(myMouse);
    glutMotionFunc(myMouseMove);

    myInit(lookat);

    // Creat Menu for Car View (top,bottom,right,left)
    // Car View
    glutCreateMenu(CarViewMenu);
    glutAddMenuEntry("larger", 6);
    glutAddMenuEntry("smaller", 7);
    glutAddMenuEntry("top", 1);
    glutAddMenuEntry("bottom", 2);
    glutAddMenuEntry("right", 3);
    glutAddMenuEntry("left", 4);
    glutAddMenuEntry("normal", 5);


    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glutMainLoop();
}


