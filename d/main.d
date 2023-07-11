/// smallpt, a Path Tracer by Kevin Beason, originally written in 2008.
/// Ported to D.
module main;

import std.algorithm;
import std.math;
import std.file;
import std.stdio;
import std.conv;

import smallpt.primitives;
import smallpt.vec;
import smallpt.rng;

static immutable Sphere[9] spheres = [
    // Left
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(0.75, .25, .25), ReflectionType.Diffuse),
    // Right
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, 0.75), ReflectionType.Diffuse),
    // Back
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(0.75, 0.75, 0.75), ReflectionType.Diffuse),
    // Front
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), ReflectionType.Diffuse),
    // Bottom
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(0.75, 0.75, 0.75), ReflectionType.Diffuse),
    // Top
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(0.75, 0.75, 0.75), ReflectionType.Diffuse),
    // Mirror
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * 0.999, ReflectionType.Specular),
    // Glass
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * 0.999, ReflectionType.Reflective),
    // Light
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), ReflectionType.Diffuse),
];

Vec radiance(in Ray ray, const int depth, scope ref Rng rng) @trusted @nogc nothrow {
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

        const Vec direction = () {
            Vec result = u * cos(r1) * r2s;

            result += nl.cross(u) * sin(r1) * r2s;

            result += nl * sqrt(1 - r2);

            return result.normalize();
        }();

        return sphere.emission + f * radiance(Ray(x, direction), depth + 1, rng);
    case ReflectionType.Specular:
        // Ideal SPECULAR reflection
        const Vec rad = radiance(Ray(x, ray.direction - (n * 2 * n.dot(ray.direction))), depth + 1, rng);
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
            return sphere.emission + (f * radiance(reflectionRay, depth + 1, rng));
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

        const Vec rad = () {
            if (depth > 2) {
                // Russian roulette
                if (rng.advanceDouble() < P) {
                    return radiance(reflectionRay, depth + 1, rng) * RP;
                }

                return radiance(Ray(x, tdir), depth + 1, rng) * TP;
            }

            return radiance(reflectionRay, depth + 1, rng) * Re +
                radiance(Ray(x, tdir), depth + 1, rng) * Tr;
        }();

        return sphere.emission + f * rad;
    }
}

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

void renderRow(
    scope ref Rng rng,
    const int y,
    scope Vec[] canvas,
    const int samples,
    const Vec cx,
    const Vec cy,
    const Ray cam,
) @trusted @nogc nothrow {
    // Loop cols
    foreach (x; 0 .. width) {
        // 2x2 subpixel rows
        foreach (sy; 0 .. 2) {
            const size_t index = (height - y - 1) * width + x;

            // 2x2 subpixel cols
            foreach (sx; 0 .. 2) {
                Vec r;

                foreach (_; 0 .. samples) {
                    const double r1 = 2 * rng.advanceDouble();
                    const double r2 = 2 * rng.advanceDouble();

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

                    // Camera rays are pushed forward to start in interior
                    const Vec position = cam.position + direction * 140;
                    r += radiance(Ray(position, direction), 0, rng) * (1.0 / samples);
                }

                canvas[index] += Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
            }
        }
    }
}

enum int width = 1024;
enum int height = 768;

static assert(height % 4 == 0);

import core.thread.osthread;

void main(string[] args) {
    const int samples = (args.length >= 2) ? to!int(args[1]) / 4 : 1;

    const cam = Ray(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize());

    const cx = Vec(width * 0.5135 / height);
    const Vec cy = cx.cross(cam.direction).normalize() * 0.5135;

    auto canvas = new Vec[](width * height);

    enum ulong randomMagic = 13889610318698649526;
    Rng rng = {randomMagic};

    // Loop over image rows
    for (int y = 0; y < height; y += 8) {
        stderr.writef!"\rRendering (%s spp) %5.2f%%"(samples * 4, 100.0 * y / (height - 1));

        Thread[8] threads = [
            new Thread(() => renderRow(rng, y, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(rng, y + 1, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(rng, y + 2, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(rng, y + 3, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(rng, y + 4, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(rng, y + 5, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(rng, y + 6, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(rng, y + 7, canvas, samples, cx, cy, cam)).start(),
        ];

        foreach (thread; threads) {
            thread.join();
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
