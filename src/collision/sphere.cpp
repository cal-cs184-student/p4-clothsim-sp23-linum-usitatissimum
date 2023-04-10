#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
  // TODO (Part 3): Handle collisions with spheres.
	Vector3D po(origin - pm.position);
	double po_d(po.norm());
	if (po_d < radius) {
		Vector3D tan_p((po_d - radius) / po_d * po + pm.position);
		Vector3D corr_vec((tan_p - pm.last_position) * (1.0 - friction));
		pm.position = pm.last_position + corr_vec;
	}
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}
