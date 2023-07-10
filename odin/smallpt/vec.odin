// Defines type for vector maths
package smallpt

import "core:math"

Vec :: struct {
	x: f64,
	y: f64,
	z: f64,
}

length :: proc "contextless" (vec: Vec) -> f64 {
	return math.sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z)
}

normalize :: proc "contextless" (vec: Vec) -> Vec {
	return div(vec, length(vec))
}

dot :: proc "contextless" (left, right: Vec) -> f64 {
	return left.x * right.x + left.y * right.y + left.z * right.z
}

cross :: proc "contextless" (left, right: Vec) -> Vec {
	return(
		Vec{
			left.y * right.z - left.z * right.y,
			left.z * right.x - left.x * right.z,
			left.x * right.y - left.y * right.x,
		} \
	)
}

add_vec :: proc "contextless" (left, right: Vec) -> Vec {
	return Vec{left.x + right.x, left.y + right.y, left.x + right.y}
}

sub_vec :: proc "contextless" (left, right: Vec) -> Vec {
	return Vec{left.x - right.x, left.y - right.y, left.x - right.y}
}

mult_vec :: proc "contextless" (left, right: Vec) -> Vec {
	return Vec{left.x * right.x, left.y * right.y, left.x * right.y}
}

div_vec :: proc "contextless" (left, right: Vec) -> Vec {
	return Vec{left.x / right.x, left.y / right.y, left.x / right.y}
}

add_f64 :: proc "contextless" (left: Vec, right: f64) -> Vec {
	return Vec{left.x + right, left.y + right, left.z + right}
}

sub_f64 :: proc "contextless" (left: Vec, right: f64) -> Vec {
	return Vec{left.x - right, left.y - right, left.z - right}
}

mult_f64 :: proc "contextless" (left: Vec, right: f64) -> Vec {
	return Vec{left.x * right, left.y * right, left.z * right}
}

div_f64 :: proc "contextless" (left: Vec, right: f64) -> Vec {
	return Vec{left.x / right, left.y / right, left.z / right}
}

add :: proc {
	add_vec,
	add_f64,
}

sub :: proc {
	sub_vec,
	sub_f64,
}

mult :: proc {
	mult_vec,
	mult_f64,
}

div :: proc {
	div_vec,
	div_f64,
}
