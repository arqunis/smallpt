module smallpt.primitives;

import smallpt.vec;
import smallpt.constants;

import std.math;

/// Material types.
enum ReflectionType {
    Diffuse,
    Specular,
    Reflective,
}

struct Sphere {
    double radius = 0.0;

    Vec position;
    Vec emission;
    Vec colour;

    ReflectionType reflection;
}

struct Ray {
    Vec position;
    Vec direction;
}

/// Calculates the distance to the intersection.
/// Returns: The distance to the ray, 0 if no hit
double intersection(in Sphere sphere, in Ray ray) @safe @nogc pure nothrow {
    // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    const Vec op = sphere.position - ray.position;

    const double b = op.dot(ray.direction);
    double det = b * b - op.dot(op) + sphere.radius * sphere.radius;
    if (det < 0) {
        return 0;
    }

    det = sqrt(det);

    double t = b - det;
    if (t > epsilon) {
        return t;
    }

    t = b + det;
    if (t > epsilon) {
        return t;
    }

    return 0;
}

bool intersection(
    scope const Sphere[] spheres,
    in Ray ray,
    scope ref double distance,
    out size_t id,
) @safe @nogc pure nothrow {
    distance = infinity;

    foreach_reverse (i, sphere; spheres) {
        double dist = intersection(sphere, ray);

        if (0 < dist && dist < distance) {
            distance = dist;
            id = i;
        }
    }

    return distance < infinity;
}
