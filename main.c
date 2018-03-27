
/*      PURPOSE : Renders any number of 3D, ray traced spheres with shading and shadowing on top of a plane.
        AUTHOR: S. Beauchemin, Tim Etchells
        DATE: 6/11/2015

        PREREQUISITES : matrix.h
 */
#define _GNU_SOURCE
#include <X11/Xlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

#define SPHERE   1
#define PLANE    2

#define EPSILON 0.00001
#define MAX_OBJECTS 100
#define MAX_INTENSITY 255.0

#define Ex 5.0 /* position of camera frame of reference */
#define Ey 5.0
#define Ez 5.0

#define Gx 0.0 /* camera gaze direction */
#define Gy 0.0
#define Gz 0.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

#define Lx 10.0 /* light position */
#define Ly -10.0
#define Lz 10.0

#define Near 1.0
#define Far 25.0

#define THETA 45.0
#define ASPECT 1.5

#define H 800 /* window height in pixels */

#define DEBUG 1

typedef struct {
    int width, height ;
} window_t ;

typedef struct {
    dmatrix_t UP ;
    dmatrix_t E ;
    dmatrix_t G ;
    dmatrix_t u, v, n ;
} camera_t ;

typedef struct {
    double r, g, b ;
} color_t ;

typedef struct {
    int type ;
    double (*intersection_function)(dmatrix_t *,dmatrix_t *) ;
    dmatrix_t M, Minv ;
    color_t specular_color, diffuse_color, ambient_color ;
    double reflectivity, specular_coeff, diffuse_coeff, ambient_coeff, f ;
} object_t ;

typedef struct {
    dmatrix_t position ;
    color_t color ;
    color_t intensity ;
} light_t ;

object_t objects[MAX_OBJECTS] ;
int nobjects = 0 ;


Display *InitX(Display *d, Window *w, int *s, window_t *Window) {

    d = XOpenDisplay(NULL) ;
    if(d == NULL) {
        printf("Cannot open display\n") ;
        exit(1) ;
    }
    *s = DefaultScreen(d) ;
    *w = XCreateSimpleWindow(d,RootWindow(d,*s),0,0,Window->width,Window->height,1,BlackPixel(d,*s),WhitePixel(d, *s)) ;
    Atom delWindow = XInternAtom(d,"WM_DELETE_WINDOW",0) ;
    XSetWMProtocols(d,*w,&delWindow,1) ;
    XSelectInput(d,*w,ExposureMask | KeyPressMask) ;
    XMapWindow(d,*w) ;
    return(d) ;
}

void SetCurrentColorX(Display *d, GC *gc, unsigned int r, unsigned int g, unsigned int b) {

    XSetForeground(d,*gc,r << 16 | g << 8 | b) ;
}

void SetPixelX(Display *d, Window w, int s, int i, int j) {

    XDrawPoint(d,w,DefaultGC(d,s),i,j) ;
}


void QuitX(Display *d, Window w) {

    XDestroyWindow(d,w) ;
    XCloseDisplay(d) ;
}

double dot_product(dmatrix_t *a, dmatrix_t *b) {

    double s  ;

    s = a->m[1][1]*b->m[1][1] + a->m[2][1]*b->m[2][1] + a->m[3][1]*b->m[3][1] ;

    return(s) ;
}

light_t *build_light(light_t *light, dmatrix_t *position, color_t color, color_t intensity) {

    dmat_alloc(&light->position,4,1) ;
    light->position = *position ;
    light->color.r = color.r ;
    light->color.g = color.g ;
    light->color.b = color.b ;
    light->intensity.r = intensity.r ;
    light->intensity.g = intensity.g ;
    light->intensity.b = intensity.b ;
    return light ;
}

window_t *build_window(window_t *Window, int height, float aspect) {

    Window->height = height ;
    Window->width =  aspect*height ;

    return(Window) ;
}

camera_t *build_camera(camera_t *Camera, window_t *Window) {

    dmat_alloc(&Camera->E,4,1) ;

    Camera->E.m[1][1] = Ex ;
    Camera->E.m[2][1] = Ey ;
    Camera->E.m[3][1] = Ez ;
    Camera->E.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->G,4,1) ;

    Camera->G.m[1][1] = Gx ;
    Camera->G.m[2][1] = Gy ;
    Camera->G.m[3][1] = Gz ;
    Camera->G.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->n,4,1) ;
    Camera->n = *dmat_normalize(dmat_sub(&Camera->E,&Camera->G)) ;
    Camera->n.l = 3 ;

    dmat_alloc(&Camera->UP,4,1) ;

    Camera->UP.l = 3 ;

    Camera->UP.m[1][1] = UPx ;
    Camera->UP.m[2][1] = UPy ;
    Camera->UP.m[3][1] = UPz ;
    Camera->UP.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->u,4,1) ;

    Camera->u = *dmat_normalize(dcross_product(&Camera->UP,&Camera->n)) ;
    Camera->v = *dmat_normalize(dcross_product(&Camera->n,&Camera->u)) ;

    return(Camera) ;
}

//returns the t-value at which the given ray intersects a generic sphere
double sphere_intersection(dmatrix_t *e, dmatrix_t *d ) {
            //a = d.x^2 + d.y^2 + d.z^2 (squared and then rooted, so these cancel out)
    double  a = pow(d->m[1][1], 2) + pow(d->m[2][1], 2) + pow(d->m[3][1], 2),
            //b = e dot d
            b = dot_product(e, d),
            //c is similar to a, but e instead of d, and subtract 1.
            c = pow(e->m[1][1], 2) + pow(e->m[2][1], 2) + pow(e->m[3][1], 2) - 1;

    double discriminant = pow(b, 2) - a*c;
    double t = -1, t2 = -1;

    if(discriminant >= 0) {
        t = -b/a + (sqrt(discriminant))/a,
        t2 = -b/a - (sqrt(discriminant))/a;
    }

    if(t2 < t)
        t = t2;

    return t;
}

//returns the t-value at which the given ray intersects a generic plane
double plane_intersection(dmatrix_t *e, dmatrix_t *d) {
    //doesn't intersect
    if(d->m[3][1] == 0)
        return -1;
    //does intersect
    else
        return -(e->m[3][1] / d->m[3][1]);
}

//returns a vector representing the intersection point, given a t value for an intersection of an object with the ray
dmatrix_t *intersection_coordinates(dmatrix_t *e, dmatrix_t *direction, double t) {
    dmatrix_t *intersection ;

    intersection = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(intersection,4,1) ;

    intersection->m[1][1] = e->m[1][1] + direction->m[1][1]*t ;
    intersection->m[2][1] = e->m[2][1] + direction->m[2][1]*t ;
    intersection->m[3][1] = e->m[3][1] + direction->m[3][1]*t ;
    intersection->m[4][1] = 1.0 ;

    return intersection ;
}

//returns the smallest t value for an i ntersection with an object
int find_min_hit_time(double t[MAX_OBJECTS]) {
    //index of smallest t-value
    int min = -1;
    for(int i = 0; i < nobjects; i++) {
        //if t[i] is nonnegative and it's a smaller value, it's the new minimum
        if(t[i] > 0 && (min == -1 || t[i] < t[min])) {
            min = i;
        }
    }
    return min;
}

dmatrix_t *ray_direction(camera_t *Camera, window_t *Window, double height, double width, int i, int j) {

    dmatrix_t *d ;

    d = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(d,4,1) ;

    d->m[1][1] = -1.0*Near*Camera->n.m[1][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[1][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[1][1] ;

    d->m[2][1] = -1.0*Near*Camera->n.m[2][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[2][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[2][1] ;

    d->m[3][1] = -1.0*Near*Camera->n.m[3][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[3][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[3][1] ;

    d->m[4][1] = 0.0 ;

    return(d) ;
}

//returns a unit vector to light source at intersection point
dmatrix_t *vector_to_light_source(dmatrix_t *intersection, dmatrix_t *light_position) {

    dmatrix_t *s ;

    s = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(s,4,1) ;

    s = dmat_sub(intersection, light_position);
    s = dmat_normalize(s);

    return s ;
}

//returns a unit vector to origin of camera frame of reference at intersection point
dmatrix_t *vector_to_center_of_projection(dmatrix_t *intersection, dmatrix_t *e) {
    dmatrix_t *v ;

    v = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(v,4,1) ;

    v = dmat_sub(intersection, e);
    v = dmat_normalize(v);

    return v ;
}

//multiplies elements in given 3D vector by given scalar
dmatrix_t *multiply_by_scalar(dmatrix_t *m, double scalar) {
    m->m[1][1] *= scalar;
    m->m[2][1] *= scalar;
    m->m[3][1] *= scalar;
    return m;
}

//returns a unit vector in the direction of the specular reflection given a surface normal and a light source
dmatrix_t *vector_to_specular_reflection(dmatrix_t *N, dmatrix_t *S) {
    dmatrix_t *r ;

    r = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(r,4,1) ;

    dmatrix_t *neg_s = multiply_by_scalar(S, -1);

    double s_dot_n = dot_product(N, S);

    double n_mag_squared = pow(N->m[1][1], 2) + pow(N->m[2][1], 2) + pow(N->m[3][1], 2);

    dmatrix_t *left_side = multiply_by_scalar(N, 2*(s_dot_n / n_mag_squared));

    r = dmat_add(neg_s, left_side);
    r = dmat_normalize(r);

    return r ;

}

color_t color_init(double r, double g, double b) {

    color_t s ;

    s.r = r ;
    s.g = g ;
    s.b = b ;

    return s ;
}

color_t color_mult(double a, color_t c) {

    color_t s ;

    s.r = a*c.r ;
    s.g = a*c.g ;
    s.b = a*c.b ;

    return s ;
}

color_t color_add(color_t c1, color_t c2) {

    color_t s ;

    s.r = c1.r + c2.r ;
    s.g = c1.g + c2.g ;
    s.b = c1.b + c2.b ;

    return s ;
}

dmatrix_t *get_normal(object_t *object, dmatrix_t *e, dmatrix_t *d, double t) {
    dmatrix_t* n = (dmatrix_t *)malloc(sizeof(dmatrix_t));
    dmat_alloc(n,4,1) ;

    if(object->type == SPHERE)
        n = intersection_coordinates(e,d,t);
    else if(object->type == PLANE) {
        n->m[1][1] = 0;
        n->m[2][1] = 0;
        n->m[3][1] = 1;
        n->m[4][1] = 1;
    }
    return n;
}
/* main ray-tracing routine. given a ray, performs the following:

for all objects in the scene {
    transforms the ray and the camera eye position with Minv of the object
    collects t values of ray intersection with generic object
}
if one or more intersections found {
    finds the intersection with minimum t value
    transforms the light with Minv
    computes the coordinates of the intersection
    using the coordinates of the intersection {
        computes the unit surface normal vector
        computes unit vector to center of the camera
        computes unit vector to light source
        computes vector of specular reflection
    }
    computes specular, diffuse, and ambient light components
    returns the R,G,B intensity of the pixel
} */
color_t shade(light_t *light,       //light object
              object_t *objects,    //array of ALL objects
              dmatrix_t *e,         //the origin of viewing system expressed in world system
              dmatrix_t *d,         //vector representing the current ray
              color_t background,   //colour of the background
              int level,            //"closeness" to the near plane for occlusion - does higher level mean closer? seems to
              int i, int j)         //coordinates of pixel (i,j)
              {

    color_t color = background;
    //array of t-value when ray intersects object
    double t_array[MAX_OBJECTS];
    //collect t-values of intersection
    for(int i = 0; i < nobjects; i++) {
        object_t *object = &objects[i];
        //tranform the ray
        dmatrix_t *e_trans = dmat_mult(&(object->Minv), e);
        dmatrix_t *d_trans = dmat_mult(&(object->Minv), d);
        //get the intersection
        t_array[i] = object->intersection_function(e_trans, d_trans);
    }

    int t_index = find_min_hit_time(t_array);
    //if there is an intersection with an object
    if(t_index != -1) {
        double t = t_array[t_index];
        object_t *object = &objects[t_index];
        //tranform light
        dmatrix_t *l_trans = dmat_mult(&(object->Minv), &(light->position));

        //compute shadow ray (from equation in notes)
        dmatrix_t *rt0 = intersection_coordinates(e, d, t);
        dmatrix_t *left_side = dmat_sub(l_trans, rt0);
        double u = object->intersection_function(rt0, left_side);

        //1 if the current pixel is in shadow
        int shadow = 0;
        //thank you floating point error
        if(u <= 1 && u >= 0.0000001) {
            shadow = 1;
        }

        color_t diffuse, specular, ambient;

        //if not in shadow
        if(!shadow) {
            //all of n, s, v, r are normalized in the functions that give them
            dmatrix_t *n = get_normal(object, e, d, t);
            dmatrix_t *s = vector_to_light_source(n, l_trans);
            dmatrix_t *v = vector_to_center_of_projection(n, e);
            dmatrix_t *r = vector_to_specular_reflection(n, s);
            //calculate diffuse, specular reflection coefficients
            double diff = object->diffuse_coeff  * dot_product(n, s);
            double spec = object->specular_coeff * pow(dot_product(v, r), object->f);

            //set diffuse, specular rgb values
            diffuse.r = (object->diffuse_color.r * light->intensity.r*diff)*MAX_INTENSITY;
            diffuse.g = (object->diffuse_color.g * light->intensity.g*diff)*MAX_INTENSITY;
            diffuse.b = (object->diffuse_color.b * light->intensity.b*diff)*MAX_INTENSITY;

            specular.r = (object->specular_color.r * light->intensity.r*spec)*MAX_INTENSITY;
            specular.g = (object->specular_color.g * light->intensity.g*spec)*MAX_INTENSITY;
            specular.b = (object->specular_color.b * light->intensity.b*spec)*MAX_INTENSITY;
        }
        else {
            diffuse.r = 0;
            diffuse.g = 0;
            diffuse.b = 0;

            specular.r = 0;
            specular.g = 0;
            specular.b = 0;
        }

        ambient.r = (object->ambient_color.r * object->ambient_coeff)*MAX_INTENSITY;
        ambient.g = (object->ambient_color.g * object->ambient_coeff)*MAX_INTENSITY;
        ambient.b = (object->ambient_color.b * object->ambient_coeff)*MAX_INTENSITY;

        color = color_add(diffuse, specular);
        color = color_add(color, ambient);
    }

    return color ;
}

object_t *build_object(int object_type, dmatrix_t *M, color_t ambient_color, color_t diffuse_color, color_t specular_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double f, double reflectivity) {

    object_t *object = malloc(sizeof(object_t)) ;

    object->M = *M ;
    dmat_alloc(&object->Minv,4,4) ;
    object->Minv = *dmat_inverse(&object->M) ;

    object->reflectivity = reflectivity ;

    object->specular_color.r = specular_color.r ;
    object->specular_color.g = specular_color.g ;
    object->specular_color.b = specular_color.b ;
    object->specular_coeff = specular_coeff ;
    object->f = f ;

    object->diffuse_color.r = diffuse_color.r ;
    object->diffuse_color.g = diffuse_color.g ;
    object->diffuse_color.b = diffuse_color.b ;
    object->diffuse_coeff = diffuse_coeff ;

    object->ambient_color.r = ambient_color.r ;
    object->ambient_color.g = ambient_color.g ;
    object->ambient_color.b = ambient_color.b ;
    object->ambient_coeff = ambient_coeff ;

    switch (object_type) {

        case SPHERE:    object->type = SPHERE ;
                        object->intersection_function = &sphere_intersection;
                        break ;
        case PLANE:     object->type = PLANE;
                        object->intersection_function = &plane_intersection;
                        break;
    }
    nobjects++ ;
    return object ;

}

int main() {

    Display *d ;
    Window w ;
    XEvent e ;

    int i, j, s ;

    camera_t Camera ;
    window_t Window ;
    light_t light ;
    dmatrix_t M, light_position ;
    color_t pixel, background, light_intensity, light_color, ambient_color, diffuse_color, specular_color ;
    double height, width, aspect, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity ;

    /* Set the background color */

    background.r = 0.0 ;
    background.g = 0.0 ;
    background.b = 0.0 ;

    /* Set up light position, intensity, and color */

    dmat_alloc(&light_position,4,1) ;

    light_position.m[1][1] = Lx ;
    light_position.m[2][1] = Ly ;
    light_position.m[3][1] = Lz ;
    light_position.m[4][1] = 1.0 ;

    light_intensity.r = 1.0 ;
    light_intensity.g = 1.0 ;
    light_intensity.b = 1.0 ;

    light_color.r = 1.0 ;
    light_color.g = 1.0 ;
    light_color.b = 1.0 ;

    light = *build_light(&light,&light_position,light_color,light_intensity) ;

    /* Build display window and synthetic camera */

    Window = *build_window(&Window,H,ASPECT) ;
    Camera = *build_camera(&Camera,&Window) ;

    /* Create a first sphere */

    dmat_alloc(&M,4,4) ;
    M = *dmat_identity(&M) ;

    //prompt user to set up spheres
    //don't pick spec/diffuse colours, that takes too long...
    if(!DEBUG) {
        printf("Enter how many spheres you would like to render: ");
        int num_spheres;
        scanf("%d", &num_spheres);

        reflectivity = 0.2 ;

        specular_color.r = 1.0 ;
        specular_color.g = 1.0 ;
        specular_color.b = 1.0 ;
        specular_coeff = 0.4 ;
        f = 10.0 ;

        diffuse_color.r = 1.0 ;
        diffuse_color.g = 1.0 ;
        diffuse_color.b = 1.0 ;
        diffuse_coeff = 0.4 ;

        ambient_coeff = 0.2 ;

        for(int i = 0; i < num_spheres; i++) {
            double r, g, b;
            printf("Enter x position of sphere %d: ", i);
            scanf("%lf", &M.m[1][4]);
            printf("Enter y position of sphere %d: ", i);
            scanf("%lf", &M.m[2][4]);
            printf("Enter z position of sphere %d: ", i);
            scanf("%lf", &M.m[3][4]);
            printf("Enter ambient red of sphere %d: ", i);
            scanf("%lf", &r);
            printf("Enter ambient green of sphere %d: ", i);
            scanf("%lf", &g);
            printf("Enter ambient blue of sphere %d: ", i);
            scanf("%lf", &b);

            objects[i] = *build_object(SPHERE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity);
        }
    }

    else {
        M.m[1][4] = 0.0 ; /* sphere at center of the world coordinates */
        M.m[2][4] = 0.0 ;
        M.m[3][4] = 1;

        reflectivity = 0.2 ;

        specular_color.r = 1.0 ;
        specular_color.g = 1.0 ;
        specular_color.b = 1.0 ;
        specular_coeff = 0.4 ;
        f = 10.0 ;

        diffuse_color.r = 1.0 ;
        diffuse_color.g = 1.0 ;
        diffuse_color.b = 1.0 ;
        diffuse_coeff = 0.4 ;

        ambient_color.r = 1.0 ;
        ambient_color.g = 0.1 ;
        ambient_color.b = 0.1 ;
        ambient_coeff = 0.2 ;

        objects[nobjects] = *build_object(SPHERE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;
    }

    M = *dmat_identity(&M) ;
    M.m[3][4] = 0.0 ;           //plane at z = 0

    ambient_color.r = 0.9 ;
    ambient_color.g = 0.8 ;
    ambient_color.b = 0.05 ;
    ambient_coeff = 0.2 ;

    objects[nobjects] = *build_object(PLANE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;

    /* Set near plane dimensions */

    aspect = ASPECT ;
    height = Near*tan(M_PI/180.0 * THETA/2.0) ;
    width = height*aspect ;

    dmatrix_t direction ;
    dmat_alloc(&direction,4,1) ;

    d = InitX(d,&w,&s,&Window) ;
    XNextEvent(d, &e) ;

    while (1) {
        XNextEvent(d, &e) ;
        if (e.type == Expose) {

            /* ray tracer main loop */

            for (i = 0 ; i < Window.width ; i++) {
                for (j = 0  ; j < Window.height ; j++) {
                    direction = *ray_direction(&Camera,&Window,height,width,i,j) ;
                    pixel = shade(&light,objects,&Camera.E,&direction,background,2,i,j) ;
                    SetCurrentColorX(d,&(DefaultGC(d,s)),(int)pixel.r,(int)pixel.g,(int)pixel.b) ;
                    SetPixelX(d,w,s,i,Window.height - (j + 1)) ;
                }
            }
        }
        if (e.type == KeyPress) {
            //causes the program to exit when 'q' is pressed
        	char c[255];
        	KeySym key;
        	//gives a warning (implicit declaration of XLookupString), but it works
        	if(XLookupString(&e.xkey, c, 255, &key)) {
        		if(c[0] == 'q')
        			break;
        	}
        }
        if (e.type == ClientMessage)
            break ;
    }
    free_dmatrix(direction.m,1,direction.l,1,direction.c) ;
    QuitX(d,w) ;
}
