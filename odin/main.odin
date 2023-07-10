/// smallpt, a Path Tracer by Kevin Beason, originally written in 2008.
/// Ported to Odin.
package main

import "core:fmt"
import "core:strconv"
import "core:math"

import "smallpt"

foreign import libc "system:c"

foreign libc {
	erand48 :: proc "c" (Xi: [^]u16) -> f64 ---
}

spheres := [?]smallpt.Sphere {
	// Left
	smallpt.Sphere{
		1e5,
		smallpt.Vec{1e5 + 1, 40.8, 81.6},
		smallpt.Vec{},
		smallpt.Vec{0.75, .25, .25},
		smallpt.ReflectionType.Diffuse,
	},
	// Right
	smallpt.Sphere{
		1e5,
		smallpt.Vec{-1e5 + 99, 40.8, 81.6},
		smallpt.Vec{},
		smallpt.Vec{.25, .25, 0.75},
		smallpt.ReflectionType.Diffuse,
	},
	// Back
	smallpt.Sphere{
		1e5,
		smallpt.Vec{50, 40.8, 1e5},
		smallpt.Vec{},
		smallpt.Vec{0.75, 0.75, 0.75},
		smallpt.ReflectionType.Diffuse,
	},
	// Front
	smallpt.Sphere{
		1e5,
		smallpt.Vec{50, 40.8, -1e5 + 170},
		smallpt.Vec{},
		smallpt.Vec{},
		smallpt.ReflectionType.Diffuse,
	},
	// Bottom
	smallpt.Sphere{
		1e5,
		smallpt.Vec{50, 1e5, 81.6},
		smallpt.Vec{},
		smallpt.Vec{0.75, 0.75, 0.75},
		smallpt.ReflectionType.Diffuse,
	},
	// Top
	smallpt.Sphere{
		1e5,
		smallpt.Vec{50, -1e5 + 81.6, 81.6},
		smallpt.Vec{},
		smallpt.Vec{0.75, 0.75, 0.75},
		smallpt.ReflectionType.Diffuse,
	},
	// Mirror
	smallpt.Sphere{
		16.5,
		smallpt.Vec{27, 16.5, 47},
		smallpt.Vec{},
		smallpt.mult(smallpt.Vec{1, 1, 1}, 0.999),
		smallpt.ReflectionType.Specular,
	},
	// Glass
	smallpt.Sphere{
		16.5,
		smallpt.Vec{73, 16.5, 78},
		smallpt.Vec{},
		smallpt.mult(smallpt.Vec{1, 1, 1}, 0.999),
		smallpt.ReflectionType.Reflective,
	},
	// Light
	smallpt.Sphere{
		600,
		smallpt.Vec{50, 681.6 - .27, 81.6},
		smallpt.Vec{12, 12, 12},
		smallpt.Vec{},
		smallpt.ReflectionType.Diffuse,
	},
}

/*
Vec radiance(in Ray ray, const int depth, ref ushort[3] Xi) @trusted @nogc nothrow {
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

        if (erand48(Xi) >= p) {
            // Russian Roulette
            return sphere.emission;
        }

        f *= 1 / p;
    }

    switch (sphere.reflection) {
    case smallpt.ReflectionType.Diffuse:
        // Ideal DIFFUSE reflection
        const double r1 = 2 * PI * erand48(Xi);
        const double r2 = erand48(Xi);
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

        return sphere.emission + f * radiance(Ray(x, direction), depth + 1, Xi);
    case smallpt.ReflectionType.Specular:
        // Ideal SPECULAR reflection
        const Vec rad = radiance(Ray(x, ray.direction - (n * 2 * n.dot(ray.direction))), depth + 1, Xi);
        return sphere.emission + f * rad;
    case smallpt.ReflectionType.Reflective:
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
            return sphere.emission + (f * radiance(reflectionRay, depth + 1, Xi));
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
                if (erand48(Xi) < P) {
                    return radiance(reflectionRay, depth + 1, Xi) * RP;
                }

                return radiance(Ray(x, tdir), depth + 1, Xi) * TP;
            }

            return radiance(reflectionRay, depth + 1, Xi) * Re +
                radiance(Ray(x, tdir), depth + 1, Xi) * Tr;
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
    const int y,
    scope Vec[] canvas,
    const int samples,
    const Vec cx,
    const Vec cy,
    const Ray cam,
) @trusted @nogc nothrow {
    ushort[3] Xi = [0, 0, cast(ushort)(y * y * y)];

    // Loop cols
    foreach (x; 0 .. width) {
        // 2x2 subpixel rows
        foreach (sy; 0 .. 2) {
            const size_t index = (height - y - 1) * width + x;

            // 2x2 subpixel cols
            foreach (sx; 0 .. 2) {
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

                    // Camera rays are pushed forward to start in interior
                    const Vec position = cam.position + direction * 140;
                    r += radiance(Ray(position, direction), 0, Xi) * (1.0 / samples);
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

    // Loop over image rows
    for (int y = 0; y < height; y += 8) {
        stderr.writef!"\rRendering (%s spp) %5.2f%%"(samples * 4, 100.0 * y / (height - 1));

        Thread[8] threads = [
            new Thread(() => renderRow(y, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(y + 1, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(y + 2, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(y + 3, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(y + 4, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(y + 5, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(y + 6, canvas, samples, cx, cy, cam)).start(),
            new Thread(() => renderRow(y + 7, canvas, samples, cx, cy, cam)).start(),
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
*/
