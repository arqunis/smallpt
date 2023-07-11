module smallpt.radiance;

import smallpt.vec;
import smallpt.primitives;
import smallpt.rng;

import std.math : sqrt, sin, cos, fabs, PI;
import std.algorithm : max;

public Vec radiance(
    scope const Sphere[] spheres,
    in Ray ray,
    const int depth,
    scope ref Rng rng,
) @safe @nogc pure nothrow {
    const result = intersection(spheres, ray);
    if (!result.hit) {
        // If missed, return black
        return Vec();
    }

    const double distance = result.distance;
    // The hit object
    const Sphere* sphere = result.sphere;

    const Vec x = ray.position + ray.direction * distance;
    const Vec n = (x - sphere.position).normalize();
    const Vec nl = (n.dot(ray.direction) < 0) ? n : (n * -1);

    Vec f = sphere.colour;

    if (depth > 5) {
        // max reflection
        const double p = max(f.x, f.y, f.z);

        if (rng.advanceDouble() >= p) {
            // Russian Roulette
            return sphere.emission;
        }

        f *= 1 / p;
    }

    switch (sphere.reflection) {
    case ReflectionType.Diffuse:
        // Ideal DIFFUSE reflection
        const double r1 = 2 * PI * rng.advanceDouble();
        const double r2 = rng.advanceDouble();
        const double r2s = sqrt(r2);

        const Vec u = () {
            auto result = (fabs(nl.x) > 0.1) ? Vec(0, 1) : Vec(1);
            return result.cross(nl).normalize();
        }();

        const Vec direction = {
            Vec result = u * cos(r1) * r2s;

            result += nl.cross(u) * sin(r1) * r2s;

            result += nl * sqrt(1 - r2);

            return result.normalize();
        }();

        return sphere.emission + f * radiance(spheres, Ray(x, direction), depth + 1, rng);
    case ReflectionType.Specular:
        // Ideal SPECULAR reflection
        const Vec rad = radiance(spheres, Ray(x, ray.direction - (n * 2 * n.dot(ray.direction))), depth + 1, rng);
        return sphere.emission + f * rad;
    case ReflectionType.Reflective:
    default:
        // Ideal dielectric REFRACTION
        const reflectionRay = Ray(x, ray.direction - (n * 2 * n.dot(ray.direction)));

        // Ray from outside going in?
        const bool into = n.dot(nl) > 0;

        const double nnt = into ? (1 / 1.5) : 1.5;
        const double ddn = ray.direction.dot(nl);
        const double cos2t = 1 - (nnt * nnt) * (1 - ddn * ddn);

        // Total internal reflection
        if (cos2t < 0) {
            return sphere.emission + (f * radiance(spheres, reflectionRay, depth + 1, rng));
        }

        const Vec tdir = (ray.direction * nnt -
                n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))))
            .normalize();

        const double c = 1 - (into ? -ddn : tdir.dot(n));

        const double Re = 0.04 + (0.96 * c * c * c * c * c);
        const double Tr = 1 - Re;
        const double P = 0.25 + 0.05 * Re;
        const double RP = Re / P;
        const double TP = Tr / (1 - P);

        const Vec rad = {
            if (depth > 2) {
                // Russian roulette
                if (rng.advanceDouble() < P) {
                    return radiance(spheres, reflectionRay, depth + 1, rng) * RP;
                }

                return radiance(spheres, Ray(x, tdir), depth + 1, rng) * TP;
            }

            return radiance(spheres, reflectionRay, depth + 1, rng) * Re +
                radiance(spheres, Ray(x, tdir), depth + 1, rng) * Tr;
        }();

        return sphere.emission + f * rad;
    }
}
