module smallpt.rng;

import std.math;

public ulong squirrel13(const ulong seed) @safe @nogc pure nothrow {
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

public struct Rng {
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
