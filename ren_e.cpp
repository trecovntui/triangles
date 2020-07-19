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

#include <Eigen/Dense>
#include <Eigen/LU>

#include <iostream>
#include <vector>

#define CLIP_W 100
#define CLIP_W 100
#define IMG_W 300
#define IMG_H 300
#define DISTS -5

typedef Eigen::Matrix3d Tr;
typedef Eigen::Vector2d Vector_2;
typedef Eigen::Vector3d Vector_3;

Vector_3 translate(0, 0, -30);

Tr scale = (Tr() <<
	100,   0,   0,
	  0, 100,   0,
    0,   0, 100
).finished();

Tr rotate = (Tr() <<
	0.866,   0.5, 0,
	 -0.5, 0.866, 0,
	    0,     0, 1
).finished();

// invert y
Tr invert_y = (Tr() <<
		1,  0, 0,
		0, -1, 0,
		0,  0, 1
).finished();

Vector_3 cube[8] = {
	Vector_3(0, 0, 0),
	Vector_3(1, 0, 0),
	Vector_3(1, 1, 0),
	Vector_3(0, 1, 0),
	Vector_3(0, 0, 1),
	Vector_3(1, 0, 1),
	Vector_3(1, 1, 1),
	Vector_3(0, 1, 1)
};

Vector_3 triangle[3] = {
	Vector_3(0, 0, -0.03),
	Vector_3(1, 0, 0.01),
	Vector_3(0, 1, 0)
};

// 222,  69,  72
//  22, 124,  76
//  38, 128, 167

#define NUM_T 1

/*
Vector_3 triang[NUM_T * 3][2] = {
	{ Vector_3(0, 0, -60), Vector_3(222,  69,  72) },
	{ Vector_3(100, 0,  0), Vector_3(222,  69,  72) },
	{ Vector_3(0, 100,  0), Vector_3(222,  69,  72) },

	{ Vector_3(0, 0, -30), Vector_3( 38, 128, 167) },
	{ Vector_3(100, 0, -30), Vector_3( 38, 128, 167) },
	{ Vector_3(0, 100, -30), Vector_3( 38, 128, 167) }
};
*/

Vector_3 triang[NUM_T * 3][2] = {
	{ Vector_3(0, 0, 0), Vector_3(222,  69,  72) },
	{ Vector_3(200, 0,  0), Vector_3( 22, 124,  76) },
	{ Vector_3(0, 200,  0), Vector_3(38, 128, 167) }
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

  int n_faces = 0;

  // Step 1: Project Primitives
  int v_index = 0;
  Vector_2 v_t_p[3];
  Vector_3 tp;
  double z_i_t[3];

  Vector_2 t_pixel;

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
      tp = invert_y * tp;
      tp = rotate * tp;
      tp = translate + tp;

      // project
      v_t_p[v_index] = Vector_2(DISTS * tp.x() / tp.z(), DISTS * tp.y() / tp.z());

      // store inverted z for later computation
      z_i_t[v_index] = (1.0 / tp.z());
      v_index++;
    }

    // v_t_p has the projected triangle's vertices (2D)

    Vector_3 a_1 = Vector_3(v_t_p[0].x(), v_t_p[0].y(), 1);
    Vector_3 b_1 = Vector_3(v_t_p[1].x(), v_t_p[1].y(), 1);
    Vector_3 c_1 = Vector_3(v_t_p[2].x(), v_t_p[2].y(), 1);

    double det_t = (Tr() << a_1, b_1, c_1).finished().determinant();

    for(int i = 0; i < (IMG_W - 1); i++) {
      for(int j = 0; j < (IMG_H - 1); j++) {
        t_pixel = Vector_2(((float)((float)i - (IMG_W / 2))) / ratio, ((float)((float)j - (IMG_H / 2))) / ratio);

        Vector_3 tp_1 = Vector_3(t_pixel.x(), t_pixel.y(), 1);

        double alpha_d = (Tr() << tp_1, b_1, c_1).finished().determinant() / det_t;
        double beta_d = (Tr() << a_1, tp_1, c_1).finished().determinant() / det_t;
        double gamma_d = (Tr() << a_1, b_1, tp_1).finished().determinant() / det_t;

        if((alpha_d >= 0) && (beta_d >= 0) && (gamma_d >= 0)) {
          // pixel inside triangle

          double depth = (z_i_t[0] * alpha_d) +
                         (z_i_t[1] * beta_d) +
                         (z_i_t[2] * gamma_d);
          depth = (1.0 / depth);
          depth = (-1.0) * depth; // now depth is positive

          if(depth > zarray[i][j]) {
            Vector_3 v1c = triang[(3 * p) + 0][1];
            Vector_3 v2c = triang[(3 * p) + 1][1];
            Vector_3 v3c = triang[(3 * p) + 2][1];

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
