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

double image[IMG_W][IMG_H] = {0};

int main(int argc, char* argv[])
{
  //std::vector<Point_3> points;
  //Point_3 p;
  
  Surface_mesh sm;
  //CGAL::convex_hull_3(points.begin(), points.end(), sm);
  //CGAL::convex_hull_3(cube, cube + 8, sm);
  CGAL::convex_hull_3(triangle, triangle + 3, sm);

  //std::cout << "The convex hull contains " << num_vertices(sm) << " vertices" << std::endl;
  //std::cout << "n_faces: " << sm.number_of_faces() << std::endl;

  int n_faces = 0;

  /*
  for (Surface_mesh::Face_index face_index : sm.faces()) {
    n_faces++;
    Surface_mesh::Halfedge_index hf = sm.halfedge(face_index);
    std::cout << "face " << n_faces << ":" << std::endl;
    for(Surface_mesh::Halfedge_index hi : halfedges_around_face(hf, sm)) {
      Surface_mesh::Vertex_index vi = target(hi, sm);
      std::cout << "vi: " << sm.point(vi) << std::endl;
    }
  }

  std::cout << "n_faces: " << n_faces << std::endl;
  */

  // Step 1: Project Primitives
  int v_index = 0;
  Point_2 v_t_p[3];
  Point_3 tp;
  double z_i_t[3];

  Point_2 a, b, c; // of triangle
  Point_2 t_pixel;

  float ratio = (IMG_W / CLIP_W);

  for (Surface_mesh::Face_index face_index : sm.faces()) {
    Surface_mesh::Halfedge_index hf = sm.halfedge(face_index);
    n_faces++;
    //std::cout << "face " << n_faces << ":" << std::endl;
    v_index = 0;
    for(Surface_mesh::Halfedge_index hi : halfedges_around_face(hf, sm)) {
      Surface_mesh::Vertex_index vi = target(hi, sm);
      //std::cout << "vi: " << sm.point(vi) << std::endl;
      tp = sm.point(vi);
      tp = scale(tp);
      //tp = rotate(tp);
      tp = translate(tp);
      //std::cout << "vt: " << tp << std::endl;
      v_t_p[v_index] = Point_2(DISTS * tp.x() / tp.z(), DISTS * tp.y() / tp.z());
      z_i_t[v_index] = (1.0 / tp.z());
      int x_map = (int)((ratio * v_t_p[v_index].x()) + (IMG_W / 2));
      int y_map = (int)((ratio * v_t_p[v_index].y()) + (IMG_H / 2));
      if(((x_map < IMG_W) && (x_map >= 0)) && ((y_map < IMG_H) && (y_map >= 0))) {
	//image[x_map][y_map] = 2;
      }
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

    Vector_2 ba = Vector_2(a, b); // b - a
    Vector_2 ca = Vector_2(a, c); // c - a

    Vector_2 bc = Vector_2(c, b); // b - c
    Vector_2 ac = Vector_2(c, a); // a - c

    Vector_2 cap = ca.perpendicular(CGAL::CLOCKWISE);
    Vector_2 bap = ba.perpendicular(CGAL::CLOCKWISE);
    Vector_2 bcp = bc.perpendicular(CGAL::CLOCKWISE);

    for(int i = 0; i < (IMG_W - 1); i++) {
      for(int j = 0; j < (IMG_H - 1); j++) {
	image[i][j] = 1.0 / z_i_t[1];
        t_pixel = Point_2(((float)((float)i - (IMG_W / 2))) / ratio, ((float)((float)j - (IMG_H / 2))) / ratio);

	Vector_2 av = Vector_2(Point_2(0, 0), a);
	Vector_2 cv = Vector_2(Point_2(0, 0), c);
	Vector_2 tpv = Vector_2(Point_2(0, 0), t_pixel);
	Vector_2 pdiff_a = tpv - (2 * av);
	Vector_2 pdiff_c = tpv - (2 * cv);

	float alpha = (pdiff_a * cap) / (ba * cap);
	float beta = (pdiff_a * bap) / (ca * bap);
	float delta = (pdiff_c * bcp) / (ac * bcp);

	if((alpha >= 0) && (beta >= 0) && (delta >= 0)) {
	  // pixel inside triangle

          Vector_3 tp_1 = Vector_3(t_pixel.x(), t_pixel.y(), 1);

	  double alpha_d = CGAL::determinant(tp_1, b_1, c_1) / det_t;
	  double beta_d = CGAL::determinant(a_1, tp_1, c_1) / det_t;
	  double gamma_d = CGAL::determinant(a_1, b_1, tp_1) / det_t;

	  //std::cout << z_i_t[0] << " " << z_i_t[1] << " " << z_i_t[2] << std::endl;

	  double depth = (z_i_t[0] * alpha_d) +
	  	         (z_i_t[1] * beta_d) +
		         (z_i_t[2] * gamma_d);
	  //std::cout << "inv_depth: " << depth << std::endl;
	  depth = (1.0 / depth);
	  //std::cout << "depth: " << depth << std::endl;

	  image[i][j] = depth;
	}
      }
    }
    // pixel loop done

  }

  for(int i = 0; i < (IMG_W - 1); i++) {
    for(int j = 0; j < (IMG_H - 1); j++) {
      std::cout << image[i][j] << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
