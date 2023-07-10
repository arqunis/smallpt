/// smallpt, a Path Tracer by Kevin Beason, originally written in 2008.
/// Ported to D.
module main;

import std.algorithm;
import std.math;
import std.file;
import std.stdio;
import std.conv;

import core.sys.posix.stdlib : erand48;

import smallpt.constants;
import smallpt.primitives;
import smallpt.vec;

static immutable Sphere[9] spheres = [
    // Left
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), ReflectionType.Diffuse),
    // Right
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), ReflectionType.Diffuse),
    // Back
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), ReflectionType.Diffuse),
    // Front
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), ReflectionType.Diffuse),
    // Bottom
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), ReflectionType.Diffuse),
    // Top
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), ReflectionType.Diffuse),
    // Mirror
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, ReflectionType.Specular),
    // Glass
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, ReflectionType.Reflective),
    // Light
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), ReflectionType.Diffuse),
];

Vec radiance(in Ray ray, int depth, ref ushort[3] Xi) @trusted @nogc nothrow {
    // Distance to intersection
    double distance = infinity;
    // Id of intersected object
    size_t id = 0;

    if (!intersection(spheres, ray, distance, id)) {
        // If missed, return black
        return Vec();
    }

    // The hit object
    const Sphere* sphere = &spheres[id];

    Vec x = ray.position + ray.direction * distance;
    Vec n = (x - sphere.position).normalize();
    Vec nl = n.dot(ray.direction) < 0 ? n : n * -1;
    Vec f = sphere.colour;

    // max reflection
    const double p = max(f.x, f.y, f.z);

    if (++depth > 5) {
        if (erand48(Xi) >= p) {
            // R.R.
            return sphere.emission;
        }

        f = f * (1 / p);
    }

    switch (sphere.reflection) {
        // Ideal DIFFUSE reflection
    case ReflectionType.Diffuse: {
            const double r1 = 2 * PI * erand48(Xi);
            const double r2 = erand48(Xi);
            const double r2s = sqrt(r2);

            Vec w = nl;
            Vec u = () {
                if (fabs(w.x) > 0.1) {
                    return Vec(0, 1).cross(w).normalize();
                }

                return Vec(1).cross(w).normalize();
            }();
            Vec v = w.cross(u);
            Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s +
                    w * sqrt(1 - r2))
                .normalize();

            return sphere.emission + f * radiance(Ray(x, d), depth, Xi);
        }
        // Ideal SPECULAR reflection
    case ReflectionType.Specular: {
            Vec rad = radiance(
                Ray(x, ray.direction - n * 2 * n.dot(ray.direction)), depth, Xi);
            return sphere.emission + f * rad;
        }
        // Ideal dielectric REFRACTION
    case ReflectionType.Reflective:
    default: {
            auto reflection_ray = Ray(x, ray.direction - n * 2 * n.dot(ray.direction));

            // Ray from outside going in?
            bool into = n.dot(nl) > 0;

            double nc = 1;
            double nt = 1.5;
            double nnt = into ? (nc / nt) : (nt / nc);
            double ddn = ray.direction.dot(nl);
            double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);

            // Total internal reflection
            if (cos2t < 0) {
                return sphere.emission + f * radiance(reflection_ray, depth, Xi);
            }

            Vec tdir = (ray.direction * nnt -
                    n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))))
                .normalize();

            double a = nt - nc;
            double b = nt + nc;
            double R0 = a * a / b * b;
            double c = 1 - (into ? -ddn : tdir.dot(n));

            double Re = R0 + (1 - R0) * c * c * c * c * c;
            double Tr = 1 - Re;
            double P = 0.25 + 0.05 * Re;
            double RP = Re / P;
            double TP = Tr / (1 - P);

            Vec rad = () {
                if (depth > 2) {
                    // Russian roulette
                    if (erand48(Xi) < P) {
                        return radiance(reflection_ray, depth, Xi) * RP;
                    }

                    return radiance(Ray(x, tdir), depth, Xi) * TP;
                }

                return radiance(reflection_ray, depth, Xi) * Re +
                    radiance(Ray(x, tdir), depth, Xi) * Tr;
            }();

            return sphere.emission + f * rad;
        }
    }
}

enum int width = 1024;
enum int height = 768;

T clamp(T)(in T x) @safe @nogc pure nothrow {
    if (x < 0) {
        return 0;
    }

    if (x > 1) {
        return 1;
    }

    return x;
}

int toInt(const double x) @safe @nogc pure nothrow {
    return cast(int)(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

void main(string[] args) {
    const int samples = (args.length >= 2) ? to!int(args[1]) / 4 : 1;

    const auto cam = Ray(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize());

    const auto cx = Vec(width * 0.5135 / height);
    const Vec cy = cx.cross(cam.direction).normalize() * 0.5135;

    Vec[] canvas = new Vec[](width * height);

    // Loop over image rows
    foreach (y; 0 .. height) {
        stderr.writef!"\rRendering (%s spp) %5.2f%%"(samples * 4, 100. * y / (height - 1));

        // Loop cols
        ushort[3] Xi = [0, 0, cast(ushort)(y * y * y)];

        foreach (x; 0 .. width) {
            // 2x2 subpixel rows
            for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++) {
                // 2x2 subpixel cols
                for (int sx = 0; sx < 2; sx++) {
                    Vec r;

                    foreach (_; 0 .. samples) {
                        const double r1 = 2 * erand48(Xi);
                        const double r2 = 2 * erand48(Xi);

                        const double dx =
                            r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        const double dy =
                            r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

                        const Vec direction = () {
                            Vec result;
                            result =
                                cx * (((sx + 0.5 + dx) / 2 + x) / width - 0.5);
                            result +=
                                cy * (((sy + 0.5 + dy) / 2 + y) / height - 0.5);
                            result += cam.direction;
                            return result.normalize();
                        }();

                        // Camera rays are pushed forward to start in
                        // interior
                        const Vec position = cam.position + direction * 140;
                        r += radiance(Ray(position, direction), 0, Xi) * (1.0 / samples);
                    }

                    canvas[i] += Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
                }
            }
        }
    }

    stderr.writeln();
    stderr.flush();

    // Write image to PPM file.
    auto f = File("image.ppm", "wb");

    f.writefln!"P3\n%s %s\n%s\n"(width, height, 255);

    foreach (i; 0 .. (width * height)) {
        f.writef!"%s %s %s "(toInt(canvas[i].x), toInt(canvas[i].y), toInt(canvas[i].z));
    }

    f.flush();
}
