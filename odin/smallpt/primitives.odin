package smallpt

import "core:math"

// Material types.
ReflectionType :: enum u8 {
    Diffuse,
    Specular,
    Reflective,
}

Sphere :: struct {
    radius: f64,

    position: Vec,
    emission: Vec,
    colour: Vec,

    reflection: ReflectionType,
}

Ray :: struct {
    position: Vec,
    direction: Vec,
}

// Calculates the distance to the intersection.
// Returns the distance to the ray, nil if no hit
@(private="file")
intersection_sphere :: proc(sphere: Sphere, ray: Ray) -> Maybe(f64) {
    // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    op := sub(sphere.position , ray.position)

    b := dot(op,ray.direction)
    det := (b * b) - dot(op, op) + (sphere.radius * sphere.radius);
    if det < 0 {
        return 0;
    }

    det = math.sqrt(det);

    t := b - det;
    if t > epsilon {
        return t;
    }

    t = b + det;
    if t > epsilon {
        return t;
    }

    return nil;
}

IntersectionResult :: struct {
    distance: f64,
    sphere: ^Sphere,
}

intersection :: proc(spheres: []Sphere, ray: Ray) -> (IntersectionResult, bool) {
    distance := infinity;
    id := 0;

    #reverse for _, i in spheres {
        #partial switch v in intersection_sphere(spheres[i], ray) {
        case f64:
            if v < distance {
                distance = v;
                id = i;
            }
        }
    }

    return IntersectionResult{distance, &spheres[id]}, distance < infinity;
}
