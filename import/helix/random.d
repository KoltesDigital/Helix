/*
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

    Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

    Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

    Neither name of Victor Nakoryakov nor the names of
    its contributors may be used to endorse or promote products
    derived from this software without specific prior written
    permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.

Copyright (C) 2006. Victor Nakoryakov.

OpenMesh/D maintainance by Bill Baxter licensed under the LGPL v2.1.
Maintenance by Jonathan Giroux (Bloutiouf) licensed under the LGPL v2.1.
http://www.gnu.org/licenses/lgpl-2.1.html
*/
/**
Module with random number generators.

Authors:
	Victor Nakoryakov (nail-mail[at]mail.ru),
	Bill Baxter,
	Jonathan Giroux (Bloutiouf)
*/
module helix.random;

import helix.basic;

version (unittest)
{
	/**
	Repeat the expression several times.
	Hopefully is the expression always true...
	*/
	private void loop(lazy void exp)
	{
		foreach (i; 0..100)
			exp();
	}
}

struct Rand48Engine
{
	private
	{
		static immutable
		{
			ulong a = 25214903917;
			ulong b = 1L;
			ulong m = 1uL << 48;
			ulong mask = m - 1;
		}
		
		ulong x = 0;
	}
	
	static immutable
	{
		uint min = 0;
		uint max = uint.max;
	}
	
	/**
	Generates next pseudo-random number.
	Returns:
		Pseudo-random number in closed interval [this.min; this.max]
	*/
	pure nothrow @safe uint pop()
	{
		x = (a * x + b) & mask;
		return cast(uint)(x >> 16);
	}
	
	/**
	Reinitializes engine. Sets new _seed used for pseudo-random number generation.
	
	If two different linear congruential engines are initialized with the same
	_seed they will generate equivalent pseudo-number sequences.
	Params:
		nx = New _seed used for pseudo-random numbers generation.
	*/
	pure nothrow @safe @property void seed(ulong nx)
	{
		x = nx & mask;
	}
	
	unittest
	{
		Rand48Engine e1;
		e1.seed = 12345;
		loop(e1.pop());
		
		Rand48Engine e2 = e1;
		
		// must generate the same sequence
		loop(assert(e1.pop() == e2.pop()));
		
		e1.seed = 54321;
		e2.seed = 54321;
		
		// must generate the same sequence
		loop(assert(e1.pop() == e2.pop()));
	}
}

/*********************************************************************/
struct MersenneTwisterEngine
{
	private
	{
		uint[n] s = void;
		uint next = 0;
		
		static pure nothrow @safe uint twist(uint m, uint u, uint v)
		{
			return m ^ (((u & 0x8000_0000u) | (v & 0x7fff_ffffu)) >> 1) ^
				(-(u & 0x1u) & 0x9908_b0dfu);
		}
	}
	
	static immutable
	{
		uint min = 0;
		uint max = uint.max;
	
		uint n = 624;
		uint m = 397;
	}
	
	pure nothrow @trusted uint pop()
	{
		if (next >= n) // overflow, engine reload needed
		{
			uint* p = s.ptr;
			
			for (int i = n - m; i--; ++p)
				*p = twist( p[m], p[0], p[1] );
			
			for (int i = m; --i; ++p)
				*p = twist( p[m - n], p[0], p[1] );
			
			*p = twist( p[m - n], p[0], s[0] );
			
			next = 0;
		}
		
		// use 'state.ptr[next]' instead of 'state[next]' to
		// suppress array bound checks, namely performance penalty
		uint x = s.ptr[next];
		++next;
		
		x ^= (x >> 11);
		x ^= (x <<  7) & 0x9d2c_5680u;
		x ^= (x << 15) & 0xefc6_0000u;
		return x ^ (x >> 18);
	}
	
	pure nothrow @safe @property void seed(uint x)
	{
		s[0] = x;
		for (int i = 1; i < n; ++i)
			s[i] = 1_812_433_253u * (s[i-1] ^ (s[i-1] >> 30)) + i;
		
		next = 1;
	}
	
	unittest
	{
		MersenneTwisterEngine e1;
		e1.seed = 12345;
		loop(e1.pop());
		
		MersenneTwisterEngine e2 = e1;
		
		// must generate the same sequence
		loop(assert(e1.pop() == e2.pop()));
		
		e1.seed = 54321;
		e2.seed = 54321;
		
		// must generate the same sequence
		loop(assert(e1.pop() == e2.pop()));
	}
}

/********************************************************************/
struct UnitUniformEngine(BaseEngine, bool closedLeft, bool closedRight)
{
	private
	{
		static immutable
		{
			real range = cast(real)(BaseEngine.max - BaseEngine.min);
			real increment = (BaseEngine.max > uint.max) ? 2L : 0.2L;
			real denominator = range + (closedLeft ? 0 : increment) + (closedRight ? 0 : increment);
		}
		
		BaseEngine baseEngine;
	}
	
	static immutable
	{
		real min = (closedLeft ? 0 : increment) * (1/denominator);
		real max = (range + (closedLeft ? 0 : increment)) * (1/denominator);
	}
	
	pure nothrow @safe real pop()
	{
		auto x = baseEngine.pop();
		
		static if (
			is (typeof(BaseEngine.pop) : real) && // base engine pops float-type values
			cast(real)BaseEngine.min == this.min &&
			cast(real)BaseEngine.max == this.max)
		{
			// no manipulations required, return value as is.
			return cast(real)x;
		}
		else
		{
			return (cast(real)(x - BaseEngine.min) + (closedLeft ? 0 : increment)) * (1/denominator);
		}
	}
	
	pure nothrow @safe @property void seed(uint x)
	{
		baseEngine.seed = x;
	}
	
	unittest
	{
		alias UnitUniformEngine!(Rand48Engine, true, true) fullClosed;
		alias UnitUniformEngine!(Rand48Engine, true, false) closedLeft;
		alias UnitUniformEngine!(Rand48Engine, false, true) closedRight;
		alias UnitUniformEngine!(Rand48Engine, false, false) fullOpened;
		
		static assert(fullClosed.min == 0.0L);
		static assert(fullClosed.max == 1.0L);
		
		static assert(closedLeft.min == 0.0L);
		static assert(closedLeft.max < 1.0L);
		
		static assert(closedRight.min > 0.0L);
		static assert(closedRight.max == 1.0L);
		
		static assert(fullOpened.min > 0.0L);
		static assert(fullOpened.max < 1.0L);
	}
}

/********************************************************************/
struct HighresUnitUniformEngine(BaseEngine, bool closedLeft, bool closedRight)
{
	private
	{
		static immutable 
		{
			real rawMax = uint.max * 0x1p32 + uint.max;
			real increment = 2.0L;
			real denominator = rawMax + (closedLeft ? 0 : increment) + (closedRight ? 0 : increment);
		}
		
		BaseEngine baseEngine;
	}
	
	static immutable
	{
		real min = (closedLeft ? 0 : increment) * (1 / denominator);
		real max = (rawMax + (closedLeft ? 0 : increment)) * (1 / denominator);
	}
	
	pure nothrow @safe real pop()
	{
		static if (
			is (typeof(BaseEngine.pop) : real) && // base engine pops float-type values
			cast(real)BaseEngine.min == this.min &&
			cast(real)BaseEngine.max == this.max)
		{
			// no manipulations required, return value as is.
			return cast(real)baseEngine.pop();
		}
		else
		{
			// this is necessary condition to generate truly uniform
			// result. However it is possible to use base engine with any range,
			// but this feature isn't implemented for now and can be introduced
			// in future.
			static assert( BaseEngine.min == 0 && BaseEngine.max == uint.max );
			
			uint a = cast(uint)baseEngine.pop();
			uint b = cast(uint)baseEngine.pop();
			return (a * 0x1p32 + b + (closedLeft ? 0 : increment)) * (1 / denominator);
		}
	}
	
	pure nothrow @safe @property void seed(uint x)
	{
		baseEngine.seed = x;
	}
	
	unittest
	{
		alias HighresUnitUniformEngine!(Rand48Engine, true, true) fullClosed;
		alias HighresUnitUniformEngine!(Rand48Engine, true, false) closedLeft;
		alias HighresUnitUniformEngine!(Rand48Engine, false, true) closedRight;
		alias HighresUnitUniformEngine!(Rand48Engine, false, false) fullOpened;
		
		static assert(fullClosed.min == 0.0L);
		static assert(fullClosed.max == 1.0L);
		
		static assert(closedLeft.min == 0.0L);
		static assert(closedLeft.max < 1.0L);
		
		static assert(closedRight.min > 0.0L);
		static assert(closedRight.max == 1.0L);
		
		static assert(fullOpened.min > 0.0L);
		static assert(fullOpened.max < 1.0L);
	}
}
