module smallpt.primitives;

import smallpt.vec;

import std.math;

/// Material types.
public enum ReflectionType {
    Diffuse,
    Specular,
    Reflective,
}

public struct Sphere {
    public double radius = 0.0;

    public Vec position;
    public Vec emission;
    public Vec colour;

    public ReflectionType reflection;
}

public struct Ray {
    public Vec position;
    public Vec direction;
}

enum double epsilon = 1e-4;
enum double infinity = 1e20;

/// Calculates the distance to the intersection.
/// Returns: The distance to the ray, 0 if no hit
public double intersection(in Sphere sphere, in Ray ray) @safe @nogc pure nothrow {
    // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    const Vec op = sphere.position - ray.position;

    const double b = op.dot(ray.direction);
    double det = (b * b) - op.dot(op) + (sphere.radius * sphere.radius);
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

public struct IntersectionResult {
    public double distance;
    public const(Sphere)* sphere;
    public bool hit;
}

public IntersectionResult intersection(return const Sphere[] spheres, in Ray ray) @safe @nogc pure nothrow {
    double distance = infinity;
    size_t id = 0;

    foreach_reverse (i, sphere; spheres) {
        double dist = intersection(sphere, ray);

        if (0 < dist && dist < distance) {
            distance = dist;
            id = i;
        }
    }

    return IntersectionResult(distance, &spheres[id], distance < infinity);
}
