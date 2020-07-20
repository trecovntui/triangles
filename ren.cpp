/*

          Y      -Z
          |      /
          |     /
          |    /
          |   /
          |  /
          | /
          |/
          O----------------X
         /
        /
       /
      /
     /
    Z

 Camera is at the origin O looking in the direction -z
 Place surface at z = DISTS

*/

#include <SFML/Graphics.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Aff_transformation_3.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define CLIP_W 100
#define CLIP_W 100
#define IMG_W 500
#define IMG_H 500
#define DISTS -100

#define USE_OBJ 1

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Point_2 Point_2;
typedef CGAL::Vector_2<K> Vector_2;
typedef CGAL::Vector_3<K> Vector_3;

typedef CGAL::Aff_transformation_3<K> Transformation;

typedef struct {
  Point_3 v[3];
  Point_3 c[3];
} Triangle;

Transformation rotate(
		0.866,   0.5, 0,
		 -0.5, 0.866, 0,
		    0,     0, 1,
		    1
);

Transformation rotate_y(
		0.866,     0,  -0.5,
		    0,     1,     0,
		  0.5,     0, 0.866,
		    1
);

Transformation rotate_x(
		    1,     0,     0,
		    0, 0.866,  -0.5,
		    0,   0.5, 0.866,
		    1
);

// invert y
Transformation invert_y(
		1,  0, 0,
		0, -1, 0,
		0,  0, 1,
    1
);

//// Red: 222,  69,  72
// Green:  22, 124,  76
/// Blue:  38, 128, 167

Point_3 c_green( 22, 124,  76);
Point_3 c_beige(245, 245, 220);

// light
Vector_3 light(0, 0, 1);

#define NUM_T 1

/*
Point_3 triang[NUM_T * 3][2] = {
	{ Point_3(  0,   0, -60), Point_3(222,  69,  72) },
	{ Point_3(100,   0,   0), Point_3(222,  69,  72) },
	{ Point_3(  0, 100,   0), Point_3(222,  69,  72) },

	{ Point_3(  0,   0, -30), Point_3( 38, 128, 167) },
	{ Point_3(100,   0, -30), Point_3( 38, 128, 167) },
	{ Point_3(  0, 100, -30), Point_3( 38, 128, 167) }
};
*/

Point_3 triang[NUM_T * 3][2] = {
	{ Point_3(  0,   0,  0), Point_3(222,  69,  72) },
	{ Point_3(200,   0,  0), Point_3( 22, 124,  76) },
	{ Point_3(  0, 200,  0), Point_3( 38, 128, 167) }
};


sf::Uint8* pixels = new sf::Uint8[IMG_W * IMG_H * 4];
double zarray[IMG_W][IMG_H] = {0};

double determine(Vector_3 x, Vector_3 y, Vector_3 z) {
  double deter = 0;

  deter += x.x() * (y.y() - z.y());
  deter += y.x() * (z.y() - x.y());
  deter += z.x() * (x.y() - y.y());

  return deter;
}

std::vector <Triangle> triangles;

void read_obj(std::string obj_file) {
  std::vector <Point_3> vertices;

  std::fstream object_file;
  object_file.open(obj_file, std::ios::in);
  if(object_file.is_open()) {
    std::string line;
    while(getline(object_file, line)) {
      // okay
      std::stringstream line_s(line);

      if(!line_s.good()) {
        continue;
      }

      std::string what;

      line_s >> what;

      if(what.compare("v") == 0) {
        // it's a vertex
        double a, b, c;
        line_s >> a; line_s >> b; line_s >> c;
        vertices.push_back(Point_3(a, b, c));
      }

      if(what.compare("f") == 0) {
        // it's a polygon
        std::vector <int> vertis;

        while(line_s.good()) {
          std::string subs;
          line_s >> subs;

          if(!(subs.length() > 0)) {
            continue;
          }

          std::stringstream subs_s(subs);
          // just need the first
          std::string v_index_s;
          getline(subs_s, v_index_s, '/');

          vertis.push_back(std::stoi(v_index_s) - 1);
        }

        if(vertis.size() == 4) {
          // It's a quadrilateral
          // Need to break it up into triangles
          Point_3 qv[4];
          for(int i = 0; i < 4; i++) {
            qv[i] = vertices[vertis[i]];
          }
          Vector_3 s1 = Vector_3(qv[0], qv[1]);
          Vector_3 s2 = Vector_3(qv[0], qv[2]);
          Vector_3 s3 = Vector_3(qv[0], qv[3]);

          if((CGAL::cross_product(s1, s2) * CGAL::cross_product(s1, s3)) < 0) {
            // s1 is diagonal
            Triangle t1 = {
              { qv[0], qv[1], qv[2] },
              { c_green, c_green, c_green }
            };
            Triangle t2 = {
              { qv[0], qv[1], qv[3] },
              { c_green, c_green, c_green }
            };
            triangles.push_back(t1);
            triangles.push_back(t2);
          }
          else if((CGAL::cross_product(s2, s1) * CGAL::cross_product(s2, s3)) < 0) {
            // s2 is diagonal
            Triangle t1 = {
              { qv[0], qv[1], qv[2] },
              { c_green, c_green, c_green }
            };
            Triangle t2 = {
              { qv[0], qv[2], qv[3] },
              { c_green, c_green, c_green }
            };
            triangles.push_back(t1);
            triangles.push_back(t2);
          }
          else {
            // s3 is diagonal
            Triangle t1 = {
              { qv[0], qv[1], qv[3] },
              { c_green, c_green, c_green }
            };
            Triangle t2 = {
              { qv[0], qv[2], qv[3] },
              { c_green, c_green, c_green }
            };
            triangles.push_back(t1);
            triangles.push_back(t2);
          }
        } // if quadrilateral
        else if(vertis.size() == 3) {
          // triangle
          Triangle tr = {
            { vertices[vertis[0]], vertices[vertis[1]], vertices[vertis[2]] },
            { c_green, c_green, c_green }
          };
          triangles.push_back(tr);
        }
      } // what = f

    } // line while
  } // if open successful
}

int main(int argc, char* argv[])
{
  // normalise light direction
  light = light / std::sqrt(light.squared_length());

  std::cout << "n_total_triangles: " << triangles.size() << std::endl;

  sf::RenderWindow window(sf::VideoMode(IMG_W, IMG_H), "triangles", sf::Style::Titlebar | sf::Style::Close);

  Point_2 v_t_p[3];
  Point_3 tp;
  double z_i_t[3];
  Point_2 t_pixel;

  float ratio = (IMG_W / CLIP_W);

  sf::Texture texture;
  sf::Sprite sprite;
  texture.create(IMG_W, IMG_H);

  Point_3 bgcolour = c_beige;
  for(int i = 0; i < IMG_W; i++) {
    for(int j = 0; j < IMG_H; j++) {
      pixels[(j * IMG_W * 4) + (i * 4) + 0] = bgcolour.x();
      pixels[(j * IMG_W * 4) + (i * 4) + 1] = bgcolour.y();
      pixels[(j * IMG_W * 4) + (i * 4) + 2] = bgcolour.z();
      pixels[(j * IMG_W * 4) + (i * 4) + 3] = 255;
    }
  }

#if USE_OBJ

  read_obj(argv[3]);

  Transformation translate(CGAL::TRANSLATION, Vector_3(0, 0, std::stoi(argv[1])));
  Transformation scale(CGAL::SCALING, std::stoi(argv[2]));

  for(int p = 0; p < triangles.size(); p++) {
    for(int q = 0; q < 3; q++) {
      tp = triangles[p].v[q];
      // transform the object in any way
      // apply the T to tp
      tp = invert_y(tp);
      //tp = rotate_x(tp);
      tp = scale(tp);
      tp = translate(tp);

      // project
      v_t_p[q] = Point_2(DISTS * tp.x() / tp.z(), DISTS * tp.y() / tp.z());

      // store inverted z for later computation
      z_i_t[q] = (1.0 / tp.z());
    }

    Point_3 tcolour[3] = {
      triangles[p].c[0],
      triangles[p].c[1],
      triangles[p].c[2]
    };

    Vector_3 d1(triangles[p].v[0], triangles[p].v[1]);
    Vector_3 d2(triangles[p].v[0], triangles[p].v[2]);

#else // not using .obj

  Transformation translate(CGAL::TRANSLATION, Vector_3(0, 0, std::stoi(argv[1])));
  Transformation scale(CGAL::SCALING, std::stoi(argv[2]));

  for(int p = 0; p < NUM_T; p++) {
    int v_index = 0;
    for(int q = (3 * p); q < (3 * (p + 1)); q++) {
      tp = triang[q][0];
      // transform the object in any way
      // apply the T to tp
      tp = invert_y(tp);
      //tp = rotate(tp);
      tp = scale(tp);
      tp = translate(tp);

      // project
      v_t_p[v_index] = Point_2(DISTS * tp.x() / tp.z(), DISTS * tp.y() / tp.z());

      // store inverted z for later computation
      z_i_t[v_index] = (1.0 / tp.z());
      v_index++;
    }

    Point_3 tcolour[3] = {
      triang[(3 * p) + 0][1],
      triang[(3 * p) + 1][1],
      triang[(3 * p) + 2][1]
    };

    Vector_3 d1(triang[(3 * p) + 0][0], triang[(3 * p) + 1][0]);
    Vector_3 d2(triang[(3 * p) + 0][0], triang[(3 * p) + 2][0]);
#endif

    Vector_3 normal = CGAL::cross_product(d1, d2);
    normal = normal / std::sqrt(normal.squared_length());

    double dot_factor = -1 * normal * light;
    dot_factor = (dot_factor < 0) ? 0 : dot_factor;

    // v_t_p has the projected triangle's vertices (2D)

    Vector_3 a_1 = Vector_3(v_t_p[0].x(), v_t_p[0].y(), 1);
    Vector_3 b_1 = Vector_3(v_t_p[1].x(), v_t_p[1].y(), 1);
    Vector_3 c_1 = Vector_3(v_t_p[2].x(), v_t_p[2].y(), 1);

    double det_t = CGAL::determinant(a_1, b_1, c_1);

    // find bounding box
    double minx = std::min(v_t_p[0].x(), std::min(v_t_p[1].x(), v_t_p[2].x()));
    double maxx = std::max(v_t_p[0].x(), std::max(v_t_p[1].x(), v_t_p[2].x()));
    double miny = std::min(v_t_p[0].y(), std::min(v_t_p[1].y(), v_t_p[2].y()));
    double maxy = std::max(v_t_p[0].y(), std::max(v_t_p[1].y(), v_t_p[2].y()));

    int i_min = (int)((ratio * minx) + (IMG_W / 2));
    int j_min = (int)((ratio * miny) + (IMG_H / 2));
    int i_max = (int)((ratio * maxx) + (IMG_W / 2));
    int j_max = (int)((ratio * maxy) + (IMG_H / 2));

    if(i_min < 0) i_min = 0;
    if(i_max >= IMG_W) i_max = (IMG_W - 1);
    if(j_min < 0) j_min = 0;
    if(j_max >= IMG_H) j_max = (IMG_H - 1);

    for(int i = i_min; i <= i_max; i++) {
      for(int j = j_min; j <= j_max; j++) {
        t_pixel = Point_2(((float)((float)i - (IMG_W / 2))) / ratio, ((float)((float)j - (IMG_H / 2))) / ratio);

        Vector_3 tp_1 = Vector_3(t_pixel.x(), t_pixel.y(), 1);

        double alpha_d = CGAL::determinant(tp_1, b_1, c_1) / det_t;
        double beta_d = CGAL::determinant(a_1, tp_1, c_1) / det_t;
        double gamma_d = CGAL::determinant(a_1, b_1, tp_1) / det_t;

        if((alpha_d >= 0) && (beta_d >= 0) && (gamma_d >= 0)) {
          // pixel inside triangle

          double depth = (z_i_t[0] * alpha_d) +
                         (z_i_t[1] * beta_d) +
                         (z_i_t[2] * gamma_d);
          depth = (1.0 / depth);
          depth = (-1.0) * depth; // now depth is positive

          if(depth > zarray[i][j]) {
            zarray[i][j] = depth;

            Vector_3 v1c = Vector_3(Point_3(0, 0, 0), tcolour[0]);
            Vector_3 v2c = Vector_3(Point_3(0, 0, 0), tcolour[1]);
            Vector_3 v3c = Vector_3(Point_3(0, 0, 0), tcolour[2]);

            Vector_3 pix_c = (alpha_d * v1c * z_i_t[0]) +
                             (beta_d  * v2c * z_i_t[1]) +
                             (gamma_d * v3c * z_i_t[2]);
            pix_c = -1.0 * depth * dot_factor * pix_c;

            pixels[(j * IMG_W * 4) + (i * 4) + 0] = pix_c.x();
            pixels[(j * IMG_W * 4) + (i * 4) + 1] = pix_c.y();
            pixels[(j * IMG_W * 4) + (i * 4) + 2] = pix_c.z();
            pixels[(j * IMG_W * 4) + (i * 4) + 3] = 255;
          }
        }
      }
    }
    // pixel loop done
  } // triangle loop done

  texture.update(pixels);
  sprite.setTexture(texture);

  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) {
        window.close();
      }
    }
    window.clear();
    window.draw(sprite);
    window.display();
  }

  return 0;
}
