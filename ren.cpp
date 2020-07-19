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
 Place surface at z = -5 (DISTS)

*/

#include <SFML/Graphics.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <vector>

#include <CGAL/Vector_2.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Aff_transformation_3.h>

#define CLIP_W 100
#define CLIP_W 100
#define IMG_W 200
#define IMG_H 200
#define DISTS -5

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Point_2 Point_2;
typedef CGAL::Vector_2<K> Vector_2;
typedef CGAL::Vector_3<K> Vector_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

typedef CGAL::Aff_transformation_3<K> Transformation;

Transformation translate(CGAL::TRANSLATION, Vector_3(0, 0, -30));
Transformation scale(CGAL::SCALING, 100);

Transformation rotate(
		0.866,   0.5, 0,
		 -0.5, 0.866, 0,
		    0,     0, 0,
		    1
);

// invert y
Transformation invert_y(
		1,  0, 0,
		0, -1, 0,
		0,  0, 1,
    1
);

Point_3 cube[8] = {
	Point_3(0, 0, 0),
	Point_3(1, 0, 0),
	Point_3(1, 1, 0),
	Point_3(0, 1, 0),
	Point_3(0, 0, 1),
	Point_3(1, 0, 1),
	Point_3(1, 1, 1),
	Point_3(0, 1, 1)
};

Point_3 triangle[3] = {
	Point_3(0, 0, -0.03),
	Point_3(1, 0, 0.01),
	Point_3(0, 1, 0)
};

// 222,  69,  72
//  22, 124,  76
//  38, 128, 167

#define NUM_T 1

/*
Point_3 triang[NUM_T * 3][2] = {
	{ Point_3(0, 0, -60), Point_3(222,  69,  72) },
	{ Point_3(100, 0,  0), Point_3(222,  69,  72) },
	{ Point_3(0, 100,  0), Point_3(222,  69,  72) },

	{ Point_3(0, 0, -30), Point_3( 38, 128, 167) },
	{ Point_3(100, 0, -30), Point_3( 38, 128, 167) },
	{ Point_3(0, 100, -30), Point_3( 38, 128, 167) }
};
*/

Point_3 triang[NUM_T * 3][2] = {
	{ Point_3(0, 0, 0), Point_3(222,  69,  72) },
	{ Point_3(100, 0,  0), Point_3( 22, 124,  76) },
	{ Point_3(0, 100,  0), Point_3(38, 128, 167) }
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

int main(int argc, char* argv[])
{
  sf::RenderWindow window(sf::VideoMode(IMG_W, IMG_H), "triangles", sf::Style::Titlebar | sf::Style::Close);

  Surface_mesh sm;
  CGAL::convex_hull_3(triangle, triangle + 3, sm);

  int n_faces = 0;

  // Step 1: Project Primitives
  int v_index = 0;
  Point_2 v_t_p[3];
  Point_3 tp;
  double z_i_t[3];

  Point_2 a, b, c; // of projected triangle
  Point_2 t_pixel;

  float ratio = (IMG_W / CLIP_W);

  sf::Texture texture;
  sf::Sprite sprite;
  texture.create(IMG_W, IMG_H);

  for(int p = 0; p < NUM_T; p++) {
    v_index = 0;
    for(int q = (3 * p); q < (3 * (p + 1)); q++) {
      tp = triang[q][0];
      // transform the object in any way
      // apply the T to tp
      //tp = scale(tp);
      //tp = rotate(tp);
      tp = invert_y(tp);
      tp = translate(tp);
      v_t_p[v_index] = Point_2(DISTS * tp.x() / tp.z(), DISTS * tp.y() / tp.z());
      z_i_t[v_index] = (1.0 / tp.z());
      v_index++;
    }

    // v_t_p has the projected triangle's vertices (2D)
    a = v_t_p[0];
    b = v_t_p[1];
    c = v_t_p[2];

    Vector_3 a_1 = Vector_3(a.x(), a.y(), 1);
    Vector_3 b_1 = Vector_3(b.x(), b.y(), 1);
    Vector_3 c_1 = Vector_3(c.x(), c.y(), 1);

    double det_t = CGAL::determinant(a_1, b_1, c_1);

    for(int i = 0; i < (IMG_W - 1); i++) {
      for(int j = 0; j < (IMG_H - 1); j++) {
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
            Vector_3 v1c = Vector_3(Point_3(0, 0, 0), triang[(3 * p) + 0][1]);
            Vector_3 v2c = Vector_3(Point_3(0, 0, 0), triang[(3 * p) + 1][1]);
            Vector_3 v3c = Vector_3(Point_3(0, 0, 0), triang[(3 * p) + 2][1]);

            Vector_3 pix_c = (alpha_d * v1c * z_i_t[0]) +
                             (beta_d  * v2c * z_i_t[1]) +
                             (gamma_d * v3c * z_i_t[2]);
            pix_c = -1.0 * depth * pix_c;


            pixels[(j * IMG_W * 4) + (i * 4) + 0] = pix_c.x();
            pixels[(j * IMG_W * 4) + (i * 4) + 1] = pix_c.y();
            pixels[(j * IMG_W * 4) + (i * 4) + 2] = pix_c.z();
            pixels[(j * IMG_W * 4) + (i * 4) + 3] = 255;
          }
        }
      }
    }
    // pixel loop done
  }

  texture.update(pixels);
  //texture.setSmooth(true);
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
