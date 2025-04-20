/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <tarsadiy>
 * *************************
 */

#ifdef WIN32
#include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
#include <GL/gl.h>
#include <GL/glut.h>
#elif defined(__APPLE__)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#define strcasecmp _stricmp
#endif

#include <imageIO.h>

#include "vec3.h"
#include <cstdlib>
#include <ctime>
#include <random>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char *filename = NULL;

// The different display modes.
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// While solving the homework, it is useful to make the below values smaller for debugging purposes.
// The still images that you need to submit with the homework should be at the below resolution (640x480).
// However, for your own purposes, after you have solved the homework, you can increase those values to obtain higher-resolution images.
#define WIDTH 640
#define HEIGHT 480
#define MAX_RECURSION_DEPTH 5

// The field of view of the camera, in degrees.
#define fov 60.0

// Buffer to store the image when saving it to a JPEG.
unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
  bool areaLight;
  int samples;
  vec3 areaX;
  vec3 areaY;

  // Initializer for point light kept as area light source; can be changed to point light source
  Light()
  {
    areaLight = true;
    // Change samples to see the effect of area light
    samples = 16;
    areaX = vec3(1, 0, 0);
    areaY = vec3(0, 0, 1);
  }
};

struct LightRay
{
  vec3 origin;
  vec3 direction;

  LightRay() {}
  LightRay(const vec3 &ori, const vec3 &dir)
  {
    origin = ori;
    direction = dir;
  }

  const vec3 &getOrigin() const { return origin; }
  const vec3 &getDirection() const { return direction; }
  vec3 at(double t) const { return origin + t * direction; }
};

struct Material
{
  vec3 diffuse;
  vec3 specular;
  double shininess;
};
// Adjusts the specular recursive reflective component of final color : make higher to see less recursive reflection
double rec_ref_coef = 2.0;
Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

vec3 interpolateNormals(const Triangle &triangle, double u, double v)
{
  vec3 n0(triangle.v[0].normal[0], triangle.v[0].normal[1], triangle.v[0].normal[2]);
  vec3 n1(triangle.v[1].normal[0], triangle.v[1].normal[1], triangle.v[1].normal[2]);
  vec3 n2(triangle.v[2].normal[0], triangle.v[2].normal[1], triangle.v[2].normal[2]);

  double w = 1.0 - u - v;
  vec3 normal = w * n0 + u * n1 + v * n2;
  return normalise(normal);
}

vec3 interpolateColors(const Triangle &triangle, double u, double v, bool diffuse)
{
  const double *c0 = diffuse ? triangle.v[0].color_diffuse : triangle.v[0].color_specular;
  const double *c1 = diffuse ? triangle.v[1].color_diffuse : triangle.v[1].color_specular;
  const double *c2 = diffuse ? triangle.v[2].color_diffuse : triangle.v[2].color_specular;

  double w = 1.0 - u - v;
  return vec3(
      w * c0[0] + u * c1[0] + v * c2[0],
      w * c0[1] + u * c1[1] + v * c2[1],
      w * c0[2] + u * c1[2] + v * c2[2]);
}

vec3 clampColor(const vec3 &c)
{
  return vec3(std::min(1.0, c.x), std::min(1.0, c.y), std::min(1.0, c.z));
}

double interpolateShininess(const Triangle &triangle, double u, double v)
{
  double s0 = triangle.v[0].shininess;
  double s1 = triangle.v[1].shininess;
  double s2 = triangle.v[2].shininess;

  double w = 1.0 - u - v;
  return w * s0 + u * s1 + v * s2;
}

vec3 reflect(const vec3 &v, const vec3 &n)
{
  return v - 2 * dot(v, n) * n;
}

double hitSphere(vec3 center, double radius, const LightRay &ray)
{
  vec3 oc = center - ray.getOrigin();

  double a = dot(ray.getDirection(), ray.getDirection());
  double b = dot(ray.getDirection(), oc);
  double c = dot(oc, oc) - radius * radius;
  double discriminant = b * b - a * c;

  if (discriminant < 0)
  {
    return -1.0;
  }
  else
  {
    return (b - std::sqrt(discriminant)) / a;
  }
}

bool hitTriangle(const Triangle &triangle, const LightRay &ray, double &t, double &u, double &v)
{
  const double EPSILON = 1e-8;

  vec3 v0(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
  vec3 v1(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
  vec3 v2(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

  vec3 edge1 = v1 - v0;
  vec3 edge2 = v2 - v0;

  vec3 h = cross(ray.getDirection(), edge2);
  double a = dot(edge1, h);

  // Ray is almost parallel to the triangle plane
  if (fabs(a) < EPSILON)
    return false;

  double f = 1.0 / a;
  vec3 s = ray.getOrigin() - v0;
  u = f * dot(s, h);

  // Check intersection
  if (u < 0.0 || u > 1.0)
    return false;

  vec3 q = cross(s, edge1);
  v = f * dot(ray.getDirection(), q);

  // Check intersection
  if (v < 0.0 || u + v > 1.0)
    return false;

  t = f * dot(edge2, q);

  if (t > EPSILON)
  { // ray intersection
    return true;
  }
  else
  {
    return false; // this means that there is a line intersection but not a plane-ray intersection
  }
}

vec3 computeColor(const LightRay &ray, int depth)
{
  // Recursion depth limit
  if (depth > MAX_RECURSION_DEPTH)
  {
    return vec3(0, 0, 0);
  }

  double nearest_t = std::numeric_limits<double>::infinity();
  vec3 final_color(0, 0, 0);
  bool hit_anything = false;

  vec3 intersection_point;
  vec3 normal_at_intersection;
  Material material_at_intersection;

  // Spheres
  for (int i = 0; i < num_spheres; i++)
  {
    vec3 center(spheres[i].position[0], spheres[i].position[1], spheres[i].position[2]);
    double t = hitSphere(center, spheres[i].radius, ray);
    if (t > 0.0 && t < nearest_t)
    {
      nearest_t = t;
      hit_anything = true;

      // parametric to world coordinates
      intersection_point = ray.at(t);
      normal_at_intersection = normalise(intersection_point - center);

      material_at_intersection.diffuse = vec3(spheres[i].color_diffuse[0], spheres[i].color_diffuse[1], spheres[i].color_diffuse[2]);
      material_at_intersection.specular = vec3(spheres[i].color_specular[0], spheres[i].color_specular[1], spheres[i].color_specular[2]);
      material_at_intersection.shininess = spheres[i].shininess;
    }
  }

  // Traingles
  for (int i = 0; i < num_triangles; i++)
  {
    double t, u, v;
    if (hitTriangle(triangles[i], ray, t, u, v) && t < nearest_t)
    {
      nearest_t = t;
      hit_anything = true;

      // barycentric coordinates used to interpolate normals and colors at t
      intersection_point = ray.at(t);
      normal_at_intersection = interpolateNormals(triangles[i], u, v);

      material_at_intersection.diffuse = interpolateColors(triangles[i], u, v, true);
      material_at_intersection.specular = interpolateColors(triangles[i], u, v, false);
      material_at_intersection.shininess = interpolateShininess(triangles[i], u, v);
    }
  }

  // If we hit something, we need to compute the color at the nearest intersection point
  if (hit_anything)
  {
    // Ambient Light
    vec3 ambient(ambient_light[0], ambient_light[1], ambient_light[2]);
    final_color = ambient * material_at_intersection.diffuse;

    for (int i = 0; i < num_lights; i++)
    {
      Light &light = lights[i];

      double lightContribution = 0.0;

      // Determine if its an area light or point light
      int num_samples = light.areaLight ? light.samples : 1;

      // Random number generator for stratified sampling
      static bool rand_initialized = false;
      if (!rand_initialized)
      {
        srand((unsigned int)time(NULL));
        rand_initialized = true;
      }

      for (int s = 0; s < num_samples; ++s)
      {
        vec3 light_sample_pos;

        // Implementing area light
        if (light.areaLight)
        {
          // Sampling points on the area light source; sampling in minigrids in the main square 1x1 grid
          int sqrt_samples = (int)sqrt((double)light.samples);
          int u_index = s % sqrt_samples;
          int v_index = s / sqrt_samples;
          double rand_u = (u_index + ((double)rand() / RAND_MAX)) / sqrt_samples;
          double rand_v = (v_index + ((double)rand() / RAND_MAX)) / sqrt_samples;

          // Normalize the random numbers to be in the range [-0.5, 0.5]
          rand_u -= 0.5;
          rand_v -= 0.5;

          // Sample point
          light_sample_pos = vec3(light.position[0], light.position[1], light.position[2]) + rand_u * light.areaX + rand_v * light.areaY;
        }
        else
        {
          light_sample_pos = vec3(light.position[0], light.position[1], light.position[2]);
        }

        // Unit direction vector
        vec3 L = normalise(light_sample_pos - intersection_point);

        // Shadow ray with slight offset per sample
        vec3 shadowRay_origin = intersection_point + 1e-4 * L;
        LightRay shadowRay(shadowRay_origin, L);

        // Distance to light sample
        double light_distance = (light_sample_pos - intersection_point).size();

        // Shadow check
        bool in_shadow = false;

        // Shoot shadow rays
        for (int j = 0; j < num_spheres; ++j)
        {
          vec3 center(spheres[j].position[0], spheres[j].position[1], spheres[j].position[2]);
          double t_shadow = hitSphere(center, spheres[j].radius, shadowRay);
          if (t_shadow > 0.0 && t_shadow < light_distance)
          {
            in_shadow = true;
            break;
          }
        }

        if (!in_shadow)
        {
          for (int j = 0; j < num_triangles; ++j)
          {
            double t_shadow, u_shadow, v_shadow;
            if (hitTriangle(triangles[j], shadowRay, t_shadow, u_shadow, v_shadow) &&
                t_shadow < light_distance)
            {
              in_shadow = true;
              break;
            }
          }
        }

        if (!in_shadow)
        {
          lightContribution += 1.0;
        }
      }

      // Average out the light contribution over the number of samples
      lightContribution /= num_samples;

      // If not shadow, cmpute shading
      if (lightContribution > 0.0)
      {
        vec3 light_dir = normalise(vec3(light.position[0], light.position[1], light.position[2]) - intersection_point);
        // Cap dot proudct to 0
        double LdotN = std::max(0.0, dot(light_dir, normal_at_intersection));

        vec3 R = reflect(-light_dir, normal_at_intersection);
        vec3 V = normalise(-ray.getDirection());

        // Specular
        double RdotV = std::max(0.0, dot(R, V));
        double specularComp = pow(RdotV, material_at_intersection.shininess);

        // Light color
        vec3 light_color(light.color[0], light.color[1], light.color[2]);

        // Phong
        // Note that we have element-wise multiplication here i.e. seperate for each color as we have overloaded the * operator in the vec3 header
        final_color = final_color + lightContribution * light_color * (material_at_intersection.diffuse * LdotN + material_at_intersection.specular * specularComp);
      }
    }

    // Recursive Reflection

    // Compute reflection ray
    vec3 reflect_dir = reflect(normalise(ray.getDirection()), normal_at_intersection);
    LightRay reflect_ray(intersection_point + 1e-4 * reflect_dir, reflect_dir);

    // Recursively trace the reflection ray
    vec3 reflect_color = computeColor(reflect_ray, depth + 1);

    // Blend the reflection color with the local color; again element-wise multiplication for vec3
    // Note that we have to divide the specular color by the recursion coefficient to make recursive reflection effect less; this can be changed to 1.0 to see the full effect. Also this again is
    // element-wise divison.
    final_color = (vec3(1.0, 1.0, 1.0) - (material_at_intersection.specular / rec_ref_coef)) * final_color + (material_at_intersection.specular / rec_ref_coef) * reflect_color;

    // Set color to [0,1] range
    return clampColor(final_color);
  }
  else
  {
    // Background color
    return vec3(1.0, 1.0, 1.0);
  }
}

void keyboardFunc(unsigned char key, int x, int y)
{
  switch (key)
  {
  case 27:
    exit(0);
    break;
  }
}

void draw_scene()
{
  // Camera setup
  double focal_length = 1.5;
  double viewport_height = 2.0;
  double viewport_width = viewport_height * ((double)WIDTH / (double)HEIGHT); // Aspect ratio: 640x480
  vec3 camera_center = vec3(0, 0, 0);
  vec3 viewport_u = vec3(viewport_width, 0, 0);  // x
  vec3 viewport_v = vec3(0, viewport_height, 0); // y
  vec3 delta_u = viewport_u / WIDTH;
  vec3 delta_v = viewport_v / HEIGHT;
  vec3 viewport_upper_left = camera_center - vec3(0, 0, focal_length) - viewport_u / 2 + viewport_v / 2;
  vec3 pixel00 = viewport_upper_left + 0.5 * (delta_u - delta_v);

  // Anti-aliasing
  const int samples_per_pixel = 16;

  // Random number generator setup
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  for (int y = 0; y < HEIGHT; y++)
  {
    glPointSize(2.0);
    glBegin(GL_POINTS);
    for (int x = 0; x < WIDTH; x++)
    {
      vec3 pixel_color(0, 0, 0);

      // Supersampling for anti-aliasing
      for (int s = 0; s < samples_per_pixel; ++s)
      {
        // Generate random offsets within the pixel
        double rand_u = distribution(generator);
        double rand_v = distribution(generator);

        vec3 pixel_sample = pixel00 + (x + rand_u) * delta_u - (y + rand_v) * delta_v;
        vec3 ray_direction = normalise(pixel_sample - camera_center);
        LightRay ray(camera_center, ray_direction);

        pixel_color = pixel_color + computeColor(ray, 0);
      }

      // Average the color per sample
      pixel_color = pixel_color / samples_per_pixel;

      // Color between [0,1]
      pixel_color = clampColor(pixel_color);

      unsigned char red = static_cast<unsigned char>(255.999 * pixel_color.x);
      unsigned char green = static_cast<unsigned char>(255.999 * pixel_color.y);
      unsigned char blue = static_cast<unsigned char>(255.999 * pixel_color.z);

      plot_pixel(x, y, red, green, blue);
    }
    glEnd();
    glFlush();
  }
  printf("Ray tracing completed.\n");
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  // Flipping for inveted image
  y = HEIGHT - y - 1;
  plot_pixel_display(x, y, r, g, b);
  if (mode == MODE_JPEG)
    plot_pixel_jpeg(x, y, r, g, b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in saving\n");
  else
    printf("File saved successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if (strcasecmp(expected, found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parsing error; abnormal program abortion.\n");
    exit(0);
  }
}

void parse_doubles(FILE *file, const char *check, double p[3])
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check(check, str);
  fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
  printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check("rad:", str);
  fscanf(file, "%lf", r);
  printf("rad: %f\n", *r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file, "%s", s);
  parse_check("shi:", s);
  fscanf(file, "%lf", shi);
  printf("shi: %f\n", *shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv, "r");
  if (!file)
  {
    printf("Unable to open input file %s. Program exiting.\n", argv);
    exit(0);
  }

  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file, "%i", &number_of_objects);

  printf("number of objects: %i\n", number_of_objects);

  parse_doubles(file, "amb:", ambient_light);

  for (int i = 0; i < number_of_objects; i++)
  {
    fscanf(file, "%s\n", type);
    printf("%s\n", type);
    if (strcasecmp(type, "triangle") == 0)
    {
      printf("found triangle\n");
      for (int j = 0; j < 3; j++)
      {
        parse_doubles(file, "pos:", t.v[j].position);
        parse_doubles(file, "nor:", t.v[j].normal);
        parse_doubles(file, "dif:", t.v[j].color_diffuse);
        parse_doubles(file, "spe:", t.v[j].color_specular);
        parse_shi(file, &t.v[j].shininess);
      }

      if (num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if (strcasecmp(type, "sphere") == 0)
    {
      printf("found sphere\n");

      parse_doubles(file, "pos:", s.position);
      parse_rad(file, &s.radius);
      parse_doubles(file, "dif:", s.color_diffuse);
      parse_doubles(file, "spe:", s.color_specular);
      parse_shi(file, &s.shininess);

      if (num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if (strcasecmp(type, "light") == 0)
    {
      printf("found light\n");
      parse_doubles(file, "pos:", l.position);
      parse_doubles(file, "col:", l.color);

      if (num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n", type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // Hack to make it only draw once.
  static int once = 0;
  if (!once)
  {
    draw_scene();
    if (mode == MODE_JPEG)
      save_jpg();
  }
  once = 1;
}

int main(int argc, char **argv)
{
  if ((argc < 2) || (argc > 3))
  {
    printf("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if (argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if (argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc, argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(WIDTH, HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
#ifdef __APPLE__
  // This is needed on recent Mac OS X versions to correctly display the window.
  glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
#endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutKeyboardFunc(keyboardFunc);
  glutMainLoop();
}
