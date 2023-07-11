/// smallpt, a Path Tracer by Kevin Beason, originally written in 2008.
/// Ported to D.
module main;

import std.stdio;
import std.conv;
import std.math : abs, pow, sqrt, sin, cos, fabs, PI;
import std.algorithm : max;

import core.thread.osthread;

ulong squirrel13(const ulong seed) @safe @nogc pure nothrow {
    enum ulong BIT_NOISE1 = 0xB5297A4DB5297A4D;
    enum ulong BIT_NOISE2 = 0x68E31DA468E31DA4;
    enum ulong BIT_NOISE3 = 0x1B56C4E91B56C4E9;

    ulong mangled = seed;
    mangled *= BIT_NOISE1;
    mangled ^= (mangled >> 8);
    mangled += BIT_NOISE2;
    mangled ^= (mangled << 8);
    mangled *= BIT_NOISE3;
    mangled ^= (mangled >> 8);
    return mangled;
}

struct Rng {
    ulong state;

    ulong advance() @safe @nogc pure nothrow scope {
        ulong newState = squirrel13(state);
        state = newState;
        return newState;
    }

    double advanceDouble() @safe @nogc pure nothrow scope {
        return abs(cast(double) advance() / cast(double) ulong.max);
    }
}

struct Vec {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

scope:
    double length() const @safe @nogc pure nothrow {
        return sqrt(x * x + y * y + z * z);
    }

    Vec normalize() const @safe @nogc pure nothrow {
        return this / length();
    }

    double dot(in Vec right) const @safe @nogc pure nothrow {
        return x * right.x + y * right.y + z * right.z;
    }

    Vec cross(in Vec right) const @safe @nogc pure nothrow {
        return Vec(
            y * right.z - z * right.y,
            z * right.x - x * right.z,
            x * right.y - y * right.x);
    }

    Vec opBinary(string op)(in Vec right) const @safe @nogc pure nothrow {
        return Vec(
            mixin("x" ~ op ~ "right.x"),
            mixin("y" ~ op ~ "right.y"),
            mixin("z" ~ op ~ "right.z"),
        );
    }

    Vec opBinary(string op)(const double right) const @safe @nogc pure nothrow {
        return Vec(
            mixin("x" ~ op ~ "right"),
            mixin("y" ~ op ~ "right"),
            mixin("z" ~ op ~ "right"),
        );
    }

    ref Vec opOpAssign(string op)(in Vec right) @safe @nogc pure nothrow return {
        mixin("x" ~ op ~ "= right.x;");
        mixin("y" ~ op ~ "= right.y;");
        mixin("z" ~ op ~ "= right.z;");

        return this;
    }

    ref Vec opOpAssign(string op)(const double right) @safe @nogc pure nothrow return {
        mixin("x" ~ op ~ "= right;");
        mixin("y" ~ op ~ "= right;");
        mixin("z" ~ op ~ "= right;");

        return this;
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

enum double epsilon = 1e-4;
enum double infinity = 1e20;

/// Calculates the distance to the intersection.
/// Returns: The distance to the ray, 0 if no hit
double intersection(in Sphere sphere, in Ray ray) @safe @nogc pure nothrow {
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

struct IntersectionResult {
    double distance;
    const(Sphere)* sphere;
    bool hit;
}

IntersectionResult intersection(return const Sphere[] spheres, in Ray ray) @safe @nogc pure nothrow {
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

Vec radiance(
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

enum int width = 1024;
enum int height = 768;

void renderRow(
    scope ref Rng rng,
    const int y,
    scope Vec[] canvas,
    const int samples,
    in Camera camera,
) @safe @nogc pure nothrow {
    const size_t row = (height - y - 1) * width;

    // Loop cols
    foreach (x; 0 .. width) {
        const size_t index = row + x;

        // 2x2 subpixel rows
        foreach (sy; 0 .. 2) {
            // 2x2 subpixel cols
            foreach (sx; 0 .. 2) {
                Vec r;

                foreach (_; 0 .. samples) {
                    const Vec direction = {
                        Vec result;
                        result = camera.x * (((sx + 0.5) / 2 + x) / width - 0.5);
                        result += camera.y * (((sy + 0.5) / 2 + y) / height - 0.5);
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
