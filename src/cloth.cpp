#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
  // TODO (Part 1): Build a grid of masses and springs.

    for (int j = 0; j < num_height_points; j++) {
        for (int i = 0; i < num_width_points; i++) {

            // generate point mass:
            Vector3D p(0.0);
            p[0] = ((double)i / (double)num_width_points) * width;
            p[1] = ((double)j / (double)num_height_points) * height;
            p[2] = ((orientation == HORIZONTAL) ? 1.0 : ((double)rand() / (double)RAND_MAX * 0.002 - 0.001));
            if (orientation == HORIZONTAL) swap(p[1], p[2]);
            vector<int> thisrow;
            thisrow.push_back(i);
            thisrow.push_back(j);

            bool pin((std::find(pinned.begin(), pinned.end(), thisrow) != pinned.end()));
            point_masses.push_back(PointMass(p, pin));

        }
    }

    for (int j = 0; j < num_height_points; j++) {
        for (int i = 0; i < num_width_points; i++) {
            // generate structural spring:
            int index = j * num_width_points + i;

            if (i > 0) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index - 1], STRUCTURAL));
            }
            if (j > 0) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index - num_width_points], STRUCTURAL));
            }

            // generate shearing constraints:
            if (j > 0) {
                if (i > 0) {
                    springs.push_back(Spring(&point_masses[index], &point_masses[index - num_width_points - 1], SHEARING));
                }
                if (i < num_width_points - 1) {
                    springs.push_back(Spring(&point_masses[index], &point_masses[index - num_width_points + 1], SHEARING));
                }
            }

            // generate bending constraints:
            if (i > 1) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index - 2], BENDING));
            }
            if (j > 1) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index - 2 * num_width_points], BENDING));
            }
        }
    }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;


  // TODO (Part 2): Compute total force acting on each point mass.
  Vector3D force(0.0);
  for (int i = 0; i < external_accelerations.size(); i++) {
      force += external_accelerations[i] * mass;
  }

  for (int i = 0; i < point_masses.size(); i++) {
      point_masses[i].forces = force;
  }

  double F;
  Vector3D F_s;

  for (int i = 0; i < springs.size(); i++) {
      Vector3D direction(springs[i].pm_b->position - springs[i].pm_a->position);
      direction.normalize();
      switch (springs[i].spring_type) {
      case(STRUCTURAL):
          if (!cp->enable_structural_constraints)
              goto skip;
          F = cp->ks * ((springs[i].pm_a->position - springs[i].pm_b->position).norm() - springs[i].rest_length);
          F_s = F * direction;
          springs[i].pm_a->forces += F_s;
          springs[i].pm_b->forces -= F_s;
          break;
      case(SHEARING):
          if (!cp->enable_structural_constraints)
              goto skip;
          F = cp->ks * ((springs[i].pm_a->position - springs[i].pm_b->position).norm() - springs[i].rest_length);
          F_s = F * direction;
          springs[i].pm_a->forces += F_s;
          springs[i].pm_b->forces -= F_s;
          break;
      case(BENDING):
          if (!cp->enable_structural_constraints)
              goto skip;
          F = 0.2 * cp->ks * ((springs[i].pm_a->position - springs[i].pm_b->position).norm() - springs[i].rest_length);
          F_s = F * direction;
          springs[i].pm_a->forces += F_s;
          springs[i].pm_b->forces -= F_s;
          break;
      }
  skip:
      continue;
  }

  // TODO (Part 2): Use Verlet integration to compute new point mass positions
  for (int i = 0; i < point_masses.size(); i++) {
      if (point_masses[i].pinned) continue;
      Vector3D curr_pos(point_masses[i].position);
      Vector3D last_pos(point_masses[i].last_position);
      Vector3D new_pos(curr_pos + (1.0 - cp->damping / 100.0) * (curr_pos - last_pos) + point_masses[i].forces / mass * delta_t * delta_t);

      point_masses[i].position = new_pos;
      point_masses[i].last_position = curr_pos;
  }

  // TODO (Part 4): Handle self-collisions.
  build_spatial_map();
  for (int i = 0; i < point_masses.size(); i++) {
      self_collide(point_masses[i], simulation_steps);
  }

  // TODO (Part 3): Handle collisions with other primitives.
  for (int i = 0; i < point_masses.size(); i++) {
      for (CollisionObject *co : *collision_objects) {
          co->collide(point_masses[i]);
      }
  }

  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  for (int i = 0; i < springs.size(); i++) {
      Spring s(springs[i]);
      double length((s.pm_b->position - s.pm_a->position).norm());
      Vector3D a_b(s.pm_b->position - s.pm_a->position);
      a_b.normalize();
      Vector3D b_a(-a_b);
      if (length < (1.1 * s.rest_length))
          continue;
      if ((s.pm_a->pinned) && (!(s.pm_b->pinned))) {
          s.pm_b->position += b_a * (length - s.rest_length * 1.1);
      }
      if ((s.pm_b->pinned) && (!(s.pm_a->pinned))) {
          s.pm_a->position += a_b * (length - s.rest_length * 1.1);
      }
      if ((!(s.pm_a->pinned)) && (!(s.pm_b->pinned))) {
          s.pm_b->position += b_a * (length - s.rest_length * 1.1) / 2.0;
          s.pm_a->position += a_b * (length - s.rest_length * 1.1) / 2.0;
      }
  }
}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.

  for (int i = 0; i < point_masses.size(); i++) {
      float key(hash_position(point_masses[i].position));
      if (map.find(key) == map.end()) {
          map[key] = new vector<PointMass*>;
      }
      map[key]->push_back(&point_masses[i]);
  }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.
    float key(hash_position(pm.position));
    if (map[key]->size() == 1) return;
    Vector3D buffer(0.0);
    int count(0);
    for (PointMass *candidate : *map[key]) {
        if (candidate == &pm) continue;
        Vector3D adjustment(pm.position - candidate->position);
        double dist(adjustment.norm());
        if (dist > (thickness * 2.0)) continue;
        adjustment.normalize();
        adjustment = adjustment * (thickness * 2.0 - dist);
        buffer += adjustment;
        count++;
    }
    if (count == 0) return;
    pm.position += (buffer / ((double)count) / simulation_steps);
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
    double wspacing(3.0 * width / (double)num_width_points);
    double hspacing(3.0 * height / (double )num_height_points);
    double tspacing(max(wspacing, hspacing));

    double x((pos.x - fmod(pos.x, wspacing)) / wspacing);
    double y((pos.y - fmod(pos.y, hspacing)) / hspacing);
    double z((pos.z - fmod(pos.z, tspacing)) / tspacing);

    return (x * 3.0 + y) * 3.0 + z;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
