/// smallpt, a Path Tracer by Kevin Beason, originally written in 2008.
/// Ported to D.
module main;

import smallpt.radiance;
import smallpt.primitives;
import smallpt.vec;
import smallpt.rng;

import std.math;
import std.stdio;
import std.conv;

import core.thread.osthread;

T clamp(T)(in T x) @safe @nogc pure nothrow {
    if (x < 0) {
        return 0;
    }

    if (x > 1) {
        return 1;
    }

    return x;
}

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

struct Camera {
    Ray ray;
    Vec x;
    Vec y;

    @disable this();

    this(int width, int height) @safe @nogc pure nothrow scope {
        ray = Ray(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize());

        x = Vec(width * 0.5135 / height);
        y = x.cross(ray.direction).normalize() * 0.5135;
    }
}

enum int width = 1024;
enum int height = 768;

void renderRow(
    scope ref Rng rng,
    const int y,
    scope Vec[] canvas,
    const int samples,
    in Camera camera,
) @safe @nogc pure nothrow {
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

                    const Vec direction = {
                        Vec result;
                        result =
                            camera.x * (((sx + 0.5 + dx) / 2 + x) / width - 0.5);
                        result +=
                            camera.y * (((sy + 0.5 + dy) / 2 + y) / height - 0.5);
                        result += camera.ray.direction;
                        return result.normalize();
                    }();

                    // Camera rays are pushed forward to start in interior
                    const Vec position = camera.ray.position + direction * 140;
                    r += radiance(spheres, Ray(position, direction), 0, rng) * (1.0 / samples);
                }

                canvas[index] += Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
            }
        }
    }
}

void main(string[] args) {
    const int samples = (args.length >= 2) ? to!int(args[1]) / 4 : 1;

    const camera = Camera(width, height);

    auto canvas = new Vec[](width * height);

    enum ulong randomMagic = 13889610318698649526;
    Rng rng = {randomMagic};

    static assert(height % 8 == 0);

    // Loop over image rows
    for (int y = 0; y < height; y += 8) {
        stderr.writef!"\rRendering (%s spp) %5.2f%%"(samples * 4, 100.0 * y / (height - 1));

        Thread[8] threads = [
            new Thread(() => renderRow(rng, y, canvas, samples, camera)).start(),
            new Thread(() => renderRow(rng, y + 1, canvas, samples, camera)).start(),
            new Thread(() => renderRow(rng, y + 2, canvas, samples, camera)).start(),
            new Thread(() => renderRow(rng, y + 3, canvas, samples, camera)).start(),
            new Thread(() => renderRow(rng, y + 4, canvas, samples, camera)).start(),
            new Thread(() => renderRow(rng, y + 5, canvas, samples, camera)).start(),
            new Thread(() => renderRow(rng, y + 6, canvas, samples, camera)).start(),
            new Thread(() => renderRow(rng, y + 7, canvas, samples, camera)).start(),
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
        static int toInt(const double x) @safe @nogc pure nothrow {
            return cast(int)(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
        }

        f.writef!"%s %s %s "(toInt(canvas[i].x), toInt(canvas[i].y), toInt(canvas[i].z));
    }

    f.flush();
}
