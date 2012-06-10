Helix
=====

Helix is an applied mathematics routines library for D2. It is intended to be low level-ish and useful for 3D applications developers.

Brief overview
--------------

Helix is relatively small right now, and includes following modules:

* `helix.basic`: basic routines for real numbers. Module also includes several templates that are used in other modules to provide common functions such as approximate equality, linear interpolation, etc.
* `helix.color`: includes structs for color operations in RGB (Red Green Blue), RGBA (Red Green Blue Alpha) and HSL (Hue Saturation Luminance) representations. These structs have methods that implement color-arithmetics and metods that provide ability to pack/unpack their's values to/from uint in different byte-order formats.
* `helix.linalgebra`: module mainly oriented to using in 3D graphics applications and includes 2,3,4-D vectors, 3x3 and 4x4 matrices, quaternion, and several global functions for linear algebra related manipulations.

Helix can be configured to use different real number's types used internaly by default (`float`, `double`, `real`), however you always can explicitly specify that type on variable level. So you can use `float` as internal type for first variable, `double` for second, and `real` for third.

Compilation
-----------

### Library

	dmd import/helix/*.d -oflib/helix.lib -lib -release -O -property

### Unit testing

	dmd import/helix/*.d unittesting.d -ofunittesting -O -property -unittest -g -w -wi && ./unittesting

If nothing shows up, everything is fine. Alternatively, `./unittesting` returns true (i.e. 0) when all the unit tests pass, and false otherwise, so you can use it for your automatic building systems:

	./unittesting && echo "passed" || echo "failed"

Documentation

	dmd import/helix/*.d -D -Dddoc -o-

Authors
-------

Original version by *Victor Nakoryakov* in 2006.

Maintained by *Bill Baxter* as part of the [OpenMesh/D project](http://www.dsource.org/projects/openmeshd) until 2009.

Maintained by *Jonathan Giroux (Bloutiouf)* since then.

License
-------

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.
* Neither name of Victor Nakoryakov nor the names of
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

OpenMesh/D maintainance by Bill Baxter licensed under the [LGPL v2.1](http://www.gnu.org/licenses/lgpl-2.1.html).

Maintenance by Jonathan Giroux (Bloutiouf) licensed under the [LGPL v2.1](http://www.gnu.org/licenses/lgpl-2.1.html).