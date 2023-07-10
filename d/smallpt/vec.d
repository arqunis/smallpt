/// Defines type for vector maths
module smallpt.vec;

import std.math;

public struct Vec {
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
