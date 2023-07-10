/**
 * @file main.cpp
 * @brief smallpt, a Path Tracer by Kevin Beason, originally written in 2008.
 * @details Modified to use newer features of C++ and to be more readable.
 */
#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <smallpt/constants.hpp>
#include <smallpt/primitives.hpp>
#include <smallpt/vec.hpp>

using namespace smallpt;

inline constexpr std::array<Sphere, 9> spheres = {
    Sphere(
        1e5,
        Vec(1e5 + 1, 40.8, 81.6),
        Vec(),
        Vec(.75, .25, .25),
        ReflectionType::Diffuse),  // Left
    Sphere(
        1e5,
        Vec(-1e5 + 99, 40.8, 81.6),
        Vec(),
        Vec(.25, .25, .75),
        ReflectionType::Diffuse),  // Rght
    Sphere(
        1e5,
        Vec(50, 40.8, 1e5),
        Vec(),
        Vec(.75, .75, .75),
        ReflectionType::Diffuse),  // Back
    Sphere(
        1e5,
        Vec(50, 40.8, -1e5 + 170),
        Vec(),
        Vec(),
        ReflectionType::Diffuse),  // Frnt
    Sphere(
        1e5,
        Vec(50, 1e5, 81.6),
        Vec(),
        Vec(.75, .75, .75),
        ReflectionType::Diffuse),  // Botm
    Sphere(
        1e5,
        Vec(50, -1e5 + 81.6, 81.6),
        Vec(),
        Vec(.75, .75, .75),
        ReflectionType::Diffuse),  // Top
    Sphere(
        16.5,
        Vec(27, 16.5, 47),
        Vec(),
        Vec(1, 1, 1) * .999,
        ReflectionType::Specular),  // Mirr
    Sphere(
        16.5,
        Vec(73, 16.5, 78),
        Vec(),
        Vec(1, 1, 1) * .999,
        ReflectionType::Reflective),  // Glas
    Sphere(
        600,
        Vec(50, 681.6 - .27, 81.6),
        Vec(12, 12, 12),
        Vec(),
        ReflectionType::Diffuse)  // Lite
};

static Vec
radiance(Ray const& ray, int depth, span<unsigned short> Xi) noexcept {
    // Distance to intersection
    double distance{infinity};
    // Id of intersected object
    std::size_t id{0};

    if (!intersection({spheres.data(), spheres.size()}, ray, distance, id)) {
        // If missed, return black
        return Vec();
    }

    // The hit object
    Sphere const& sphere = spheres[id];

    Vec x = ray.position + ray.direction * distance;
    Vec n = (x - sphere.position).normalize();
    Vec nl = n.dot(ray.direction) < 0 ? n : n * -1;
    Vec f = sphere.colour;

    // max reflection
    double p = std::max(std::max(f.x, f.y), f.z);

    if (++depth > 5) {
        if (erand48(Xi.data()) >= p) {
            // R.R.
            return sphere.emission;
        }

        f = f * (1 / p);
    }

    switch (sphere.reflection) {
    // Ideal DIFFUSE reflection
    case ReflectionType::Diffuse: {
        double r1 = 2 * M_PI * erand48(Xi.data());
        double r2 = erand48(Xi.data());
        double r2s = std::sqrt(r2);

        Vec w = nl;
        Vec u = [&] {
            if (std::fabs(w.x) > 0.1) {
                return Vec(0, 1).cross(w).normalize();
            }

            return Vec(1).cross(w).normalize();
        }();
        Vec v = w.cross(u);
        Vec d = (u * std::cos(r1) * r2s + v * std::sin(r1) * r2s +
                 w * std::sqrt(1 - r2))
                    .normalize();

        return sphere.emission + f * radiance(Ray(x, d), depth, Xi);
    }
    // Ideal SPECULAR reflection
    case ReflectionType::Specular: {
        Vec rad = radiance(
            Ray(x, ray.direction - n * 2 * n.dot(ray.direction)), depth, Xi);
        return sphere.emission + f * rad;
    }
    // Ideal dielectric REFRACTION
    case ReflectionType::Reflective:
    default: {
        Ray reflection_ray{x, ray.direction - n * 2 * n.dot(ray.direction)};

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
                    n * ((into ? 1 : -1) * (ddn * nnt + std::sqrt(cos2t))))
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

        Vec rad = [&] {
            if (depth > 2) {
                // Russian roulette
                if (erand48(Xi.data()) < P) {
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

inline constexpr int width = 1024;
inline constexpr int height = 768;

template <typename T>
constexpr T clamp(T x) noexcept {
    if (x < 0) {
        return 0;
    }

    if (x > 1) {
        return 1;
    }

    return x;
}

inline int to_int(double x) noexcept {
    return static_cast<int>(std::pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

int main(int argc, char* argv[]) {
    int samples = argc == 2 ? atoi(argv[1]) / 4 : 1;

    const Ray cam{Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize()};

    const Vec cx{width * 0.5135 / height};
    const Vec cy = cx.cross(cam.direction).normalize() * 0.5135;

    std::vector<Vec> canvas;
    canvas.resize(width * height);

// OpenMP
#pragma omp parallel for schedule(dynamic, 1)
    // Loop over image rows
    for (int y = 0; y < height; y++) {
        std::fprintf(
            stderr,
            "\rRendering (%d spp) %5.2f%%",
            samples * 4,
            100. * y / (height - 1));

        // Loop cols
        std::array<unsigned short, 3> Xi = {
            0, 0, static_cast<unsigned short>(y * y * y)};

        for (unsigned short x = 0; x < width; x++) {
            // 2x2 subpixel rows
            for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++) {
                // 2x2 subpixel cols
                for (int sx = 0; sx < 2; sx++) {
                    Vec r;

                    for (int s = 0; s < samples; s++) {
                        const double r1 = 2 * erand48(Xi.data());
                        const double r2 = 2 * erand48(Xi.data());

                        const double dx =
                            r1 < 1 ? std::sqrt(r1) - 1 : 1 - std::sqrt(2 - r1);
                        const double dy =
                            r2 < 1 ? std::sqrt(r2) - 1 : 1 - std::sqrt(2 - r2);

                        const Vec direction = [&] {
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
                        r += radiance(
                                 Ray(position, direction),
                                 0,
                                 {Xi.data(), Xi.size()}) *
                             (1.0 / samples);
                    }

                    canvas[i] += Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
                }
            }
        }
    }

    std::fprintf(stderr, "\n");
    std::fflush(stderr);

    // Write image to PPM file.
    std::FILE* f = std::fopen("image.ppm", "wb");

    std::fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);

    for (int i = 0; i < (width * height); i++) {
        std::fprintf(
            f,
            "%d %d %d ",
            to_int(canvas[i].x),
            to_int(canvas[i].y),
            to_int(canvas[i].z));
    }

    std::fflush(f);
    std::fclose(f);

    return 0;
}
