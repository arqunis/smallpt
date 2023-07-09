/**
 * @file main.cpp
 * @brief smallpt, a Path Tracer by Kevin Beason, originally written in 2008.
 * @details Modified to use newer features of C++ and to be more readable.
 */
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>

struct Vec {
    double x{0.0};
    double y{0.0};
    double z{0.0};

    constexpr Vec() noexcept = default;
    constexpr Vec(double x, double y, double z) noexcept : x{x}, y{y}, z{z} {}
    constexpr Vec(double x, double y) noexcept : x{x}, y{y} {}
    constexpr Vec(double x) noexcept : x{x} {}

    friend constexpr Vec operator+(Vec const& left, Vec const& right) noexcept {
        return Vec(left.x + right.x, left.y + right.y, left.z + right.z);
    }

    friend constexpr Vec operator-(Vec const& left, Vec const& right) noexcept {
        return Vec(left.x - right.x, left.y - right.y, left.z - right.z);
    }

    friend constexpr Vec operator*(Vec const& left, Vec const& right) noexcept {
        return Vec(left.x * right.x, left.y * right.y, left.z * right.z);
    }

    friend constexpr Vec operator*(Vec const& left, double right) noexcept {
        return Vec(left.x * right, left.y * right, left.z * right);
    }

    friend constexpr Vec operator/(Vec const& left, double right) noexcept {
        return Vec(left.x / right, left.y / right, left.z / right);
    }

    double length() const noexcept {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vec normalize() const noexcept {
        return *this / length();
    }

    constexpr double dot(Vec const& right) const noexcept {
        return x * right.x + y * right.y + z * right.z;
    }

    constexpr Vec cross(Vec const& right) noexcept {
        return Vec(
            y * right.z - z * right.y,
            z * right.x - x * right.z,
            x * right.y - y * right.x);
    }
};

struct Ray {
    Vec position;
    Vec direction;

    constexpr Ray(Vec position, Vec direction) noexcept
        : position{position}, direction{direction} {}
};

/// @brief Material types, which are used in `radiance`.
enum class ReflectionType { Diffuse, Specular, Reflective };

inline constexpr double eps = 1e-4;

struct Sphere {
    double radius{0.0};

    Vec position;
    Vec emission;
    Vec colour;

    ReflectionType reflection;

    constexpr Sphere(
        double radius,
        Vec position,
        Vec emission,
        Vec colour,
        ReflectionType reflection) noexcept
        : radius{radius}, position{position}, emission{emission},
          colour{colour}, reflection{reflection} {}

    /// @brief Calculates the distance to the intersection.
    /// @returns The distance to the ray, 0 if no hit
    double intersection(Ray const& ray) const noexcept {
        // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        Vec op = position - ray.position;

        double b = op.dot(ray.direction);
        double det = b * b - op.dot(op) + radius * radius;
        if (det < 0) {
            return 0;
        }

        det = std::sqrt(det);

        double t = b - det;
        if (t > eps) {
            return t;
        }

        t = b + det;
        if (t > eps) {
            return t;
        }

        return 0;
    }
};

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
    return int(std::pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

inline constexpr double inf = 1e20;

inline bool intersection(Ray const& ray, double& distance, int& id) {
    distance = inf;

    for (std::size_t i = spheres.size(); i-- > 0;) {
        double dist = spheres[i].intersection(ray);

        if (0 < dist && dist < distance) {
            distance = dist;
            id = static_cast<int>(i);
        }
    }

    return distance < inf;
}

Vec radiance(
    Ray const& ray,
    int depth,
    std::array<unsigned short, 3>& Xi) noexcept {
    // Distance to intersection
    double distance{inf};
    // Id of intersected object
    int id{0};

    if (!intersection(ray, distance, id)) {
        // If missed, return black
        return Vec();
    }

    // The hit object
    Sphere const& sphere = spheres[id];

    Vec x = ray.position + ray.direction * distance;
    Vec n = (x - sphere.position).normalize();
    Vec nl = n.dot(ray.direction) < 0 ? n : n * -1;
    Vec f = sphere.colour;

    double p = [&] {
        if (f.x > f.y && f.x > f.z) {
            return f.x;
        }

        if (f.y > f.z) {
            return f.y;
        }

        // max reflection
        return f.z;
    }();

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

int main(int argc, char* argv[]) {
    int samples = argc == 2 ? atoi(argv[1]) / 4 : 1;

    Ray cam{Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize()};

    Vec cx{width * 0.5135 / height};
    Vec cy = cx.cross(cam.direction).normalize() * 0.5135;

    Vec* canvas = new Vec[width * height];

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
                        double r1 = 2 * erand48(Xi.data());
                        double dx =
                            r1 < 1 ? std::sqrt(r1) - 1 : 1 - std::sqrt(2 - r1);

                        double r2 = 2 * erand48(Xi.data());
                        double dy =
                            r2 < 1 ? std::sqrt(r2) - 1 : 1 - std::sqrt(2 - r2);

                        Vec d =
                            cx * (((sx + 0.5 + dx) / 2 + x) / width - 0.5) +
                            cy * (((sy + 0.5 + dy) / 2 + y) / height - 0.5) +
                            cam.direction;

                        // Camera rays are pushed forward to start in
                        // interior
                        d = d.normalize();
                        Vec position = cam.position + d * 140;
                        r = r +
                            radiance(Ray(position, d), 0, Xi) * (1. / samples);
                    }

                    canvas[i] = canvas[i] +
                                Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
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

    delete[] canvas;

    return 0;
}
