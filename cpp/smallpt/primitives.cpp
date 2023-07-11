#include <smallpt/primitives.hpp>

#include <cmath>

using namespace smallpt;

static constexpr double epsilon = 1e-4;
static constexpr double infinity = 1e20;

double smallpt::intersection(Sphere const& sphere, Ray const& ray) noexcept {
    // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    const Vec op = sphere.position - ray.position;

    const double b = op.dot(ray.direction);
    double det = b * b - op.dot(op) + sphere.radius * sphere.radius;
    if (det < 0) {
        return 0;
    }

    det = std::sqrt(det);

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

bool smallpt::intersection(
    span<Sphere const> spheres,
    Ray const& ray,
    double& distance,
    std::size_t& id) noexcept {
    distance = infinity;

    for (std::size_t i = spheres.size(); i-- > 0;) {
        double dist = intersection(spheres[i], ray);

        if (0 < dist && dist < distance) {
            distance = dist;
            id = i;
        }
    }

    return distance < infinity;
}
