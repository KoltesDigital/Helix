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
Port to D2 by Jonathan Giroux (Bloutiouf) licensed under the LGPL v2.1.
http://www.gnu.org/licenses/lgpl-2.1.html
*/
/**
Module consists of basic mathematical objects oriented to working with 3D
graphics.

Those are 2,3,4-D vectors, quaternion, 3x3 and 4x4 matrices. In case of
specialization for 3D graphics there are always some features and deviation
from classical math. Here I summarize such features of helix'es linear algebra:
$(UL
    $(LI In helix paradigm of column-vector is taken. So multiplication of matrix
         by vector makes sense but multiplication of vector by matrix makes not.
         This approach conforms to rules accepted in classical math and coincides
         with OpenGL rules. However this is opposite to Direct3D paradigm where
         vector is a row. So, in helix, to combine sequence of transforms specified with
         matrices A, B, C in order A then B then C, you have to multiply them in
         back-to-front order order: M=C*B*A. )

    $(LI When an issue deal with euler angles following definitions are accepted.
         Yaw is rotation around Y, Pitch is rotaion around X, Roll is rotation around Z
         Rotations are always applied in order: Roll then Pitch then Yaw. )

    $(LI Helix matrices use column-major memory layout. I.e. matrix
        $(MAT33
            $(MAT33_ROW a, b, c)
            $(MAT33_ROW d, e, f)
            $(MAT33_ROW g, h, i)
        )
        in memory will looks like: a, d, g, b, e, h, c, f, i. This order is the same as
        in OpenGL API, but opposite to Direct3D API. However as mentioned above, Direct3D
        uses vector-row paradigm that is opposite to classical math, so D3D requires
        transposed matrix as compared to classical math to get desired transformation. As a
        result you haven't to transpose helix matrix while transmission to API even in Direct3D
        case. Normaly you haven't to remember about memory layout, just use it as in classical
        math, this feature is significant only in routines that operate with data pointers
        and plain array representation. There are reminders in such methods' documentation. )
)

Authors:
    Victor Nakoryakov (nail-mail[at]mail.ru),
	Bill Baxter,
	Jonathan Giroux (Bloutiouf)
    
Macros:
    MAT33     = <table style="border-left: double 3px #666666; border-right: double 3px #666666;
                 margin-left: 3em;">$0</table>
    MAT33_ROW = <tr><td>$1</td><td>$2</td><td>$3</td></tr>
*/
module helix.linalgebra;

import std.math;
import std.string : format;
import helix.basic,
       helix.config;

/** Defines ort names that are usualy used as indices. */
enum Ort
{
    X, ///
    Y, /// ditto
    Z, /// ditto
    W  /// ditto
}

/**
Wrapper template to provide possibility to use different float types
in implemented structs and routines.
*/
private template LinearAlgebra(float_t)
{
    private alias helix.basic.equal          equal;
    
    alias .LinearAlgebra!(float).Vector2     Vector2f;
    alias .LinearAlgebra!(float).Vector3     Vector3f;
    alias .LinearAlgebra!(float).Vector4     Vector4f;
    alias .LinearAlgebra!(float).Quaternion  Quaternionf;
    alias .LinearAlgebra!(float).Matrix22    Matrix22f;
    alias .LinearAlgebra!(float).Matrix33    Matrix33f;
    alias .LinearAlgebra!(float).Matrix44    Matrix44f;
    
    alias .LinearAlgebra!(double).Vector2    Vector2d;
    alias .LinearAlgebra!(double).Vector3    Vector3d;
    alias .LinearAlgebra!(double).Vector4    Vector4d;
    alias .LinearAlgebra!(double).Quaternion Quaterniond;
    alias .LinearAlgebra!(double).Matrix22   Matrix22d;
    alias .LinearAlgebra!(double).Matrix33   Matrix33d;
    alias .LinearAlgebra!(double).Matrix44   Matrix44d;
    
    alias .LinearAlgebra!(real).Vector2      Vector2r;
    alias .LinearAlgebra!(real).Vector3      Vector3r;
    alias .LinearAlgebra!(real).Vector4      Vector4r;
    alias .LinearAlgebra!(real).Quaternion   Quaternionr;
    alias .LinearAlgebra!(real).Matrix22     Matrix22r;
    alias .LinearAlgebra!(real).Matrix33     Matrix33r;
    alias .LinearAlgebra!(real).Matrix44     Matrix44r;

    /************************************************************************************
    Two dimensional vector.

    For formal definition of vector, meaning of possible operations and related
    information see $(LINK http://en.wikipedia.org/wiki/Vector_(spatial)).
    *************************************************************************************/
    struct Vector2
    {
        enum { length = 2u }

        align(1)
        {
            float_t x; /// Components of vector.
            float_t y; /// ditto
        }
        
        static Vector2 nan = { float_t.nan, float_t.nan }; /// Vector with both components set to NaN.
        static Vector2 zero = { 0, 0 };                    /// The zero vector 

        static Vector2 unitX = { 1, 0 };                   /// Unit vector codirectional with corresponding axis.
        static Vector2 unitY = { 0, 1 };                   /// ditto
        
        
        /**
        Method to construct vector in C-like syntax.

        Examples:
        ------------
        Vector2 myVector = Vector2(1, 2);
        ------------
        */
        pure nothrow @safe static Vector2 opCall(float_t x, float_t y)
        {
            Vector2 v;
            v.set(x, y);
            return v;
        }
        
        /** Method to construct from an array. */
        pure nothrow @safe static Vector2 opCall(float_t[2] p)
        {
            Vector2 v;
            v.set(p);
            return v;
        }
        
        /** Sets values of components to passed values. */
        pure nothrow @safe void set(float_t x, float_t y)
        {
            this.x = x;
            this.y = y;
        }
        
        /** Sets values of components to passed values. */
        pure nothrow @safe void set(float_t[2] p)
        {
            this.x = p[0];
            this.y = p[1];
        }
        
        /** Returns: Norm (also known as length, magnitude) of vector. */
        const pure nothrow @safe @property real norm()
        {
            return hypot(x, y);
        }
    
        /**
        Returns: Square of vector's norm.
        
        Since this method doesn't need calculation of square root it is better
        to use it instead of norm() when you can. For example, if you want just
        to know which of 2 vectors is longer it's better to compare their norm
        squares instead of their norm.
        */
        const pure nothrow @safe @property real normSquare()
        {
            return x*x + y*y;
        }
    
        /** Normalizes this vector. Returns: the original length. */
        pure nothrow @safe real normalize()
        {
            real len = norm();
            this /= len;
            return len;
        }
    
        /** Returns: Normalized copy of this vector. */
        const pure nothrow @safe Vector2 normalized()
        {
            real n = norm;
            return Vector2(x / n, y / n);
        }
    
        /**
        Returns: Whether this vector is unit.
        Params:
            relprec, absprec = Parameters passed to equal function while comparison of
                               norm square and 1. Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isUnit(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal( normSquare(), 1, relprec, absprec );
        }
    
        /** Returns: Axis for which projection of this vector on it will be longest. */
        const pure nothrow @safe Ort dominatingAxis()
        {
            return (x > y) ? Ort.X : Ort.Y;
        }
    
        /** Returns: Whether all components are normalized numbers. */
        const pure nothrow @safe @property bool isNormal()
        {
            return std.math.isNormal(x) && std.math.isNormal(y);
        }
    
        /** Returns: float_t pointer to x component of this vector. It's like a _ptr method for arrays. */
        pure nothrow @safe @property float_t* ptr()
        {
            return cast(float_t*)&this;
        }
    
        /** Returns: Component corresponded to passed index. */
        pure @safe float_t opIndex(size_t ort)
        in { assert(ort <= Ort.Y); }
        body
        {
            return ptr[cast(int)ort];
        }
    
        /** Assigns new _value to component corresponded to passed index. */
        pure @safe void opIndexAssign(float_t value, size_t ort)
        in { assert(ort <= Ort.Y); }
        body
        {
            ptr[cast(int)ort] = value;
        }
    
        /**
        Standard operators that have intuitive meaning, same as in classical math.
        
        Note that division operators do no cheks of value of k, so in case of division
        by 0 result vector will have infinity components. You can check this with isNormal()
        method.
        */
        const pure nothrow @safe bool opEquals(Vector2 v)
        {
            return x == v.x && y == v.y;
        }
    
        /** ditto */
        const pure nothrow @safe Vector2 opNeg()
        {
            return Vector2(-x, -y);
        }
    
        /** ditto */
        const pure nothrow @safe Vector2 opAdd(Vector2 v)
        {
            return Vector2(x + v.x, y + v.y);
        }
    
        /** ditto */
        pure nothrow @safe void opAddAssign(Vector2 v)
        {
            x += v.x;
            y += v.y;
        }
    
        /** ditto */
        const pure nothrow @safe Vector2 opSub(Vector2 v)
        {
            return Vector2(x - v.x, y - v.y);
        }
    
        /** ditto */
        pure nothrow @safe void opSubAssign(Vector2 v)
        {
            x -= v.x;
            y -= v.y;
        }
    
        /** ditto */
        const pure nothrow @safe Vector2 opMul(real k)
        {
            return Vector2(x * k, y * k);
        }
    
        /** ditto */
        pure nothrow @safe void opMulAssign(real k)
        {
            x *= k;
            y *= k;
        }
    
        /** ditto */
        const pure nothrow @safe Vector2 opMul_r(real k)
        {
            return Vector2(x * k, y * k);
        }
    
        /** ditto */
        const pure nothrow @safe Vector2 opDiv(real k)
        {
            return Vector2(x / k, y / k);
        }
    
        /** ditto */
        pure nothrow @safe void opDivAssign(real k)
        {
            x /= k;
            y /= k;
        }
    
        /** Returns: A vector perpendicular to this one */
        const pure nothrow @safe Vector2 perp() 
        {
            return Vector2(-y,x);
        }

        /** Returns: Copy of this vector with float type components */
        const pure nothrow @safe Vector2f toVector2f()
        {
            return Vector2f(cast(float)x, cast(float)y);
        }
        
        /** Returns: Copy of this vector with double type components */
        const pure nothrow @safe Vector2d toVector2d()
        {
            return Vector2d(cast(double)x, cast(double)y);
        }
        
        /** Returns: Copy of this vector with real type components */
        const pure nothrow @safe Vector2r toVector2r()
        {
            return Vector2r(cast(real)x, cast(real)y);
        }
    
        /**
        Routines known as swizzling.
        Returns:
            New vector constructed from this one and having component values
            that correspond to method name.
        */
        const pure nothrow @safe @property Vector3 xy0()    { return Vector3(x, y, 0); }
        const pure nothrow @safe @property Vector3 x0y()    { return Vector3(x, 0, y); } /// ditto

        const string toString() { return format("[",x,", ",y,"]"); }
    }
    
    /** Returns: Dot product between passed vectors. */
    pure nothrow @safe real dot(Vector2 a, Vector2 b)
    {
        return a.x * b.x + a.y * b.y;
    }
        
    /** Returns: Outer product between passed vectors. */
    pure nothrow @safe Matrix22 outer(Vector2 a, Vector2 b)
    {
        return Matrix22( a.x * b.x, a.x * b.y,
                         a.y * b.x, a.y * b.y);
    }
        
    /**
    Returns: Cross product between passed vectors. Result is scalar c
    so that [a.x a.y 0], [b.x b.y 0], [0,0,c] forms right-hand tripple.
    */
    pure nothrow @safe real cross(Vector2 a, Vector2 b)
    {
        return a.x * b.y - b.x * a.y;
    }


    alias EqualityByNorm!(Vector2).equal equal; /// Introduces approximate equality function for Vector2.
    alias Lerp!(Vector2).lerp lerp;             /// Introduces linear interpolaton function for Vector2.
    
    
    /************************************************************************************
    Three dimensional vector.

    For formal definition of vector, meaning of possible operations and related
    information see $(LINK http://en.wikipedia.org/wiki/Vector_(spatial)).
    *************************************************************************************/
    struct Vector3
    {
        enum { length = 3u }

        align(1)
        {
            float_t x; /// Components of vector.
            float_t y; /// ditto
            float_t z; /// ditto
        }
        
        static Vector3 nan = { float_t.nan, float_t.nan, float_t.nan }; /// Vector with all components set to NaN.
        static Vector3 zero = {0,0,0};    // The zero vector
        static Vector3 unitX = { 1, 0, 0 };  /// Unit vector codirectional with corresponding axis.
        static Vector3 unitY = { 0, 1, 0 };  /// ditto
        static Vector3 unitZ = { 0, 0, 1 };  /// ditto
        
        /**
        Method to construct vector in C-like syntax.

        Examples:
        ------------
        Vector3 myVector = Vector3(1, 2, 3);
        ------------
        */
        pure nothrow @safe static Vector3 opCall(float_t x, float_t y, float_t z)
        {
            Vector3 v;
            v.set(x, y, z);
            return v;
        }
        
        /** Method to construct from array */
        pure nothrow @safe static Vector3 opCall(float_t[3] p)
        {
            Vector3 v;
            v.set(p);
            return v;
        }
        
        /** Sets values of components to passed values. */
        pure nothrow @safe void set(float_t x, float_t y, float_t z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
    
        
        /** Sets values of components to passed values. */
        pure nothrow @safe void set(float_t[3] p)
        {
            this.x = p[0];
            this.y = p[1];
            this.z = p[2];
        }
    
        
        /** Returns: Norm (also known as length, magnitude) of vector. */
        const pure nothrow @safe @property real norm()
        {
            return sqrt(x*x + y*y + z*z);
        }
    
        /**
        Returns: Square of vector's norm.
        
        Since this method doesn't need calculation of square root it is better
        to use it instead of norm() when you can. For example, if you want just
        to know which of 2 vectors is longer it's better to compare their norm
        squares instead of their norm.
        */
        const pure nothrow @safe @property real normSquare()
        {
            return x*x + y*y + z*z;
        }
    
        /** Normalizes this vector. Returns: the original length */
        pure nothrow @safe real normalize()
        {
            real len = norm();
            this /= len;
            return len;
        }
    
        /** Returns: Normalized copy of this vector. */
        const pure nothrow @safe Vector3 normalized()
        {
            real n = norm;
            return Vector3(x / n, y / n, z / n);
        }
    
        /**
        Returns: Whether this vector is unit.
        Params:
            relprec, absprec = Parameters passed to equal function while comparison of
                               norm square and 1. Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isUnit(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal( normSquare(), 1, relprec, absprec );
        }
    
        /** Returns: Axis for which projection of this vector on it will be longest. */
        const pure nothrow @safe Ort dominatingAxis()
        {
            if (x > y)
                return (x > z) ? Ort.X : Ort.Z;
            else
                return (y > z) ? Ort.Y : Ort.Z;
        }
    
        /** Returns: Whether all components are normalized numbers. */
        const pure nothrow @safe bool isNormal()
        {
            return std.math.isNormal(x) && std.math.isNormal(y) && std.math.isNormal(z);
        }
    
        /** Returns: float_t pointer to x component of this vector. It's like a _ptr method for arrays. */
        pure nothrow @safe @property float_t* ptr()
        {
            return cast(float_t*)(&x);
        }
    
        /** Returns: Component corresponded to passed index. */
        pure @safe float_t opIndex(size_t ort)
        in { assert(ort <= Ort.Z); }
        body
        {
            return ptr[cast(int)ort];
        }
    
        /** Assigns new _value to component corresponded to passed index. */
        pure @safe void opIndexAssign(float_t value, size_t ort)
        in { assert(ort <= Ort.Z); }
        body
        {
            ptr[cast(int)ort] = value;
        }
    
        /**
        Standard operators that have intuitive meaning, same as in classical math.
        
        Note that division operators do no cheks of value of k, so in case of division
        by 0 result vector will have infinity components. You can check this with isNormal()
        method.
        */
        const pure nothrow @safe bool opEquals(Vector3 v)
        {
            return x == v.x && y == v.y && z == v.z;
        }
    
        /** ditto */
        const pure nothrow @safe Vector3 opNeg()
        {
            return Vector3(-x, -y, -z);
        }
    
        /** ditto */
        const pure nothrow @safe Vector3 opAdd(Vector3 v)
        {
            return Vector3(x + v.x, y + v.y, z + v.z);
        }
    
        /** ditto */
        pure nothrow @safe void opAddAssign(Vector3 v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
        }
    
        /** ditto */
        const pure nothrow @safe Vector3 opSub(Vector3 v)
        {
            return Vector3(x - v.x, y - v.y, z - v.z);
        }
    
        /** ditto */
        pure nothrow @safe void opSubAssign(Vector3 v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
        }
    
        /** ditto */
        const pure nothrow @safe Vector3 opMul(real k)
        {
            return Vector3(x * k, y * k, z * k);
        }
    
        /** ditto */
        pure nothrow @safe void opMulAssign(real k)
        {
            x *= k;
            y *= k;
            z *= k;
        }
    
        /** ditto */
        const pure nothrow @safe Vector3 opMul_r(real k)
        {
            return Vector3(x * k, y * k, z * k);
        }
    
        /** ditto */
        const pure nothrow @safe Vector3 opDiv(real k)
        {
            return Vector3(x / k, y / k, z / k);
        }
    
        /** ditto */
        pure nothrow @safe void opDivAssign(real k)
        {
            x /= k;
            y /= k;
            z /= k;
        }
        
        /** Returns: Copy of this vector with float type components */
        const pure nothrow @safe Vector3f toVector3f()
        {
            return Vector3f(cast(float)x, cast(float)y, cast(float)z);
        }
        
        /** Returns: Copy of this vector with double type components */
        const pure nothrow @safe Vector3d toVector3d()
        {
            return Vector3d(cast(double)x, cast(double)y, cast(double)z);
        }

        /** Returns: Copy of this vector with real type components */
        const pure nothrow @safe Vector3r toVector3r()
        {
            return Vector3r(cast(real)x, cast(real)y, cast(real)z);
        }

    
        /**
        Routines known as swizzling.
        Returns:
            New vector constructed from this one and having component values
            that correspond to method name.
        */
        const pure nothrow @safe @property Vector4 xyz0()        { return Vector4(x,y,z,0); }
        const pure nothrow @safe @property Vector4 xyz1()        { return Vector4(x,y,z,1); } /// ditto
        const pure nothrow @safe @property Vector2 xy()          { return Vector2(x, y); }    /// ditto
        const pure nothrow @safe @property Vector2 xz()          { return Vector2(x, z); }    /// ditto
        const pure nothrow @safe @property Vector2 yz()          { return Vector2(y, z); }    /// ditto
        
        /**
        Routines known as swizzling.
        Assigns new values to some components corresponding to method name.
        */
        pure nothrow @safe @property void xy(Vector2 v)    { x = v.x; y = v.y; }
        pure nothrow @safe @property void xz(Vector2 v)    { x = v.x; z = v.y; }        /// ditto
        pure nothrow @safe @property void yz(Vector2 v)    { y = v.x; z = v.y; }        /// ditto

        const string toString() { return format("[",x,", ",y,", ", z, "]"); }
    }
    
    /** Returns: Dot product between passed vectors. */
    pure nothrow @safe real dot(Vector3 a, Vector3 b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    
    /** Returns: Outer product between passed vectors. */
    pure nothrow @safe Matrix33 outer(Vector3 a, Vector3 b)
    {
        return Matrix33( a.x * b.x,  a.x * b.y,  a.x * b.z,
                         a.y * b.x,  a.y * b.y,  a.y * b.z,
                         a.z * b.x,  a.z * b.y,  a.z * b.z);
    }
        

    /**
    Returns: Cross product between passed vectors. Result is vector c
    so that a, b, c forms right-hand tripple.
    */
    pure nothrow @safe Vector3 cross(Vector3 a, Vector3 b)
    {
        return Vector3(
            a.y * b.z - b.y * a.z,
            a.z * b.x - b.z * a.x,
            a.x * b.y - b.x * a.y  );
    }
    
    /**
    Returns: Whether passed basis is orthogonal.
    Params:
        r, s, t =          Vectors that form a basis.
        relprec, absprec = Parameters passed to equal function while calculations.
                           Have the same meaning as in equal function.
    References:
        $(LINK http://en.wikipedia.org/wiki/Orthogonal_basis)
    */
    pure nothrow @safe bool isBasisOrthogonal(Vector3 r, Vector3 s, Vector3 t, int relprec = defrelprec, int absprec = defabsprec)
    {
        return equal( cross(r, cross(s, t)).normSquare, 0, relprec, absprec );
    }
    
    /**
    Returns: Whether passed basis is orthonormal.
    Params:
        r, s, t =   Vectors that form a basis.
        relprec, absprec = Parameters passed to equal function while calculations.
                           Have the same meaning as in equal function.
    References:
        $(LINK http://en.wikipedia.org/wiki/Orthonormal_basis)
    */
    pure nothrow @safe bool isBasisOrthonormal(Vector3 r, Vector3 s, Vector3 t, int relprec = defrelprec, int absprec = defabsprec)
    {
        return isBasisOrthogonal(r, s, t, relprec, absprec) &&
            r.isUnit(relprec, absprec) &&
            s.isUnit(relprec, absprec) &&
            t.isUnit(relprec, absprec);
    }
    
    alias EqualityByNorm!(Vector3).equal equal; /// Introduces approximate equality function for Vector3.
    alias Lerp!(Vector3).lerp lerp;             /// Introduces linear interpolation function for Vector3.
    
    /************************************************************************************
    4D vector.

    For formal definition of vector, meaning of possible operations and related
    information see $(LINK http://en.wikipedia.org/wiki/Vector_(spatial)),
    $(LINK http://en.wikipedia.org/wiki/Homogeneous_coordinates).
    *************************************************************************************/
    struct Vector4
    {
        enum { length = 4u }

        align(1)
        {
            float_t x; /// Components of vector.
            float_t y; /// ditto
            float_t z; /// ditto
            float_t w; /// ditto
        }
        
        /// Vector with all components set to NaN.
        static Vector4 nan = { float_t.nan, float_t.nan, float_t.nan, float_t.nan };
        static Vector4 zero = { 0,0,0,0 };
        static Vector4 unitX = { 1, 0, 0, 0 }; /// Unit vector codirectional with corresponding axis.
        static Vector4 unitY = { 0, 1, 0, 0 }; /// ditto
        static Vector4 unitZ = { 0, 0, 1, 0 }; /// ditto
        static Vector4 unitW = { 0, 0, 0, 1 }; /// ditto
        
        /**
        Methods to construct vector in C-like syntax.

        Examples:
        ------------
        Vector4 myVector = Vector4(1, 2, 3, 1);
        Vector4 myVector = Vector4(Vector3(1, 2, 3), 0);
        Vector4 myVector = Vector4([1,2,3,1]);
        ------------
        */
        pure nothrow @safe static Vector4 opCall(float_t x, float_t y, float_t z, float_t w)
        {
            Vector4 v;
            v.set(x, y, z, w);
            return v;
        }
        
        /** ditto */
        pure nothrow @safe static Vector4 opCall(Vector3 xyz, float_t w)
        {
            Vector4 v;
            v.set(xyz, w);
            return v;
        }
    
        /** ditto */
        pure nothrow @safe static Vector4 opCall(float_t[4] p)
        {
            Vector4 v;
            v.set(p);
            return v;
        }
    
        /** Sets values of components to passed values. */
        pure nothrow @safe void set(float_t x, float_t y, float_t z, float_t w)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
        }
    
        /** ditto */
        pure nothrow @safe void set(Vector3 xyz, float_t w)
        {
            this.x = xyz.x;
            this.y = xyz.y;
            this.z = xyz.z;
            this.w = w;
        }
    
        /** ditto */
        pure nothrow @safe void set(float_t[4] p)
        {
            this.x = p[0];
            this.y = p[1];
            this.z = p[2];
            this.w = p[3];
        }
    
        

        /**
        Returns: Norm (also known as length, magnitude) of vector.
        
        W-component is taken into account.
        */
        const pure nothrow @safe @property real norm()
        {
            return sqrt(x*x + y*y + z*z + w*w);
        }
    
        /**
        Returns: Square of vector's norm.

        W-component is taken into account.
        
        Since this method doesn't need calculation of square root it is better
        to use it instead of norm() when you can. For example, if you want just
        to know which of 2 vectors is longer it's better to compare their norm
        squares instead of their norm.
        */
        const pure nothrow @safe @property real normSquare()
        {
            return x*x + y*y + z*z + w*w;
        }
    
        /** Normalizes this vector. Returns: the original length. */
        pure nothrow @safe real normalize()
        {
            real len = norm();
            this /= len;
            return len;
        }
    
        /** Returns: Normalized copy of this vector. */
        const pure nothrow @safe Vector4 normalized()
        {
            real n = norm;
            return Vector4(x / n, y / n, z / n, w / n);
        }
    
        /**
        Returns: Whether this vector is unit.
        Params:
            relprec, absprec = Parameters passed to equal function while comparison of
                               norm square and 1. Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isUnit(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal( normSquare, 1, relprec, absprec );
        }
    
        /**
        Returns: Axis for which projection of this vector on it will be longest.
        
        W-component is taken into account.
        */
        const pure nothrow @safe Ort dominatingAxis()
        {
            if (x > y)
            {
                if (x > z)
                    return (x > w) ? Ort.X : Ort.W;
                else
                    return (z > w) ? Ort.Z : Ort.W;
            }
            else
            {
                if (y > z)
                    return (y > w) ? Ort.Y : Ort.W;
                else
                    return (z > w) ? Ort.Z : Ort.W;
            }
        }
    
        /** Returns: Whether all components are normalized numbers. */
        const pure nothrow @safe bool isNormal()
        {
            return std.math.isNormal(x) && std.math.isNormal(y) && std.math.isNormal(z) && std.math.isNormal(w);
        }
    
        /** Returns: float_t pointer to x component of this vector. It's like a _ptr method for arrays. */
        pure nothrow @safe @property float_t* ptr()
        {
            return cast(float_t*)(&x);
        }
        
        /** Returns: Component corresponded to passed index. */
        pure @safe float_t opIndex(size_t ort)
        in { assert(ort <= Ort.W); }
        body
        {
            return ptr[cast(int)ort];
        }
    
        /** Assigns new value to component corresponded to passed index. */
        pure @safe void opIndexAssign(float_t value, size_t ort)
        in { assert(ort <= Ort.W); }
        body
        {
            ptr[cast(int)ort] = value;
        }
    
        /**
        Standard operators that have intuitive meaning, same as in classical math.
        
        Note that division operators do no cheks of value of k, so in case of division
        by 0 result vector will have infinity components. You can check this with isNormal()
        method.
        */
        const pure nothrow @safe bool opEquals(Vector4 v)
        {
            return x == v.x && y == v.y && z == v.z && w == v.w;
        }
    
        /** ditto */
        const pure nothrow @safe Vector4 opNeg()
        {
            return Vector4(-x, -y, -z, -w);
        }
    
        /** ditto */
        const pure nothrow @safe Vector4 opAdd(Vector4 v)
        {
            return Vector4(x + v.x, y + v.y, z + v.z, w + v.w);
        }
    
        /** ditto */
        pure nothrow @safe void opAddAssign(Vector4 v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
            w += v.w;
        }
    
        /** ditto */
        const pure nothrow @safe Vector4 opSub(Vector4 v)
        {
            return Vector4(x - v.x, y - v.y, z - v.z, w - v.w);
        }
    
        /** ditto */
        pure nothrow @safe void opSubAssign(Vector4 v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            w -= v.w;
        }
    
        /** ditto */
        const pure nothrow @safe Vector4 opMul(real k)
        {
            return Vector4(x * k, y * k, z * k, w * k);
        }
    
        /** ditto */
        pure nothrow @safe void opMulAssign(real k)
        {
            x *= k;
            y *= k;
            z *= k;
            w *= k;
        }
    
        /** ditto */
        const pure nothrow @safe Vector4 opMul_r(real k)
        {
            return Vector4(x * k, y * k, z * k, w * k);
        }
    
        /** ditto */
        const pure nothrow @safe Vector4 opDiv(real k)
        {
            return Vector4(x / k, y / k, z / k, w / k);
        }
    
        /** ditto */
        pure nothrow @safe void opDivAssign(real k)
        {
            x /= k;
            y /= k;
            z /= k;
            w /= k;
        }
    
        /** Returns: Copy of this vector with float type components */
        const pure nothrow @safe Vector4f toVector4f()
        {
            return Vector4f(cast(float)x, cast(float)y, cast(float)z, cast(float)w);
        }
        
        /** Returns: Copy of this vector with double type components */
        const pure nothrow @safe Vector4d toVector4d()
        {
            return Vector4d(cast(double)x, cast(double)y, cast(double)z, cast(double)w);
        }
        
        /** Returns: Copy of this vector with real type components */
        const pure nothrow @safe Vector4r toVector4r()
        {
            return Vector4r(cast(real)x, cast(real)y, cast(real)z, cast(real)w);
        }
    
        /**
        Routine known as swizzling.
        Returns:
            Vector3 that has the same x, y, z components values.
        */
        const pure nothrow @safe @property Vector3 xyz()   { return Vector3(x,y,z); }    
        
        /**
        Routine known as swizzling.
        Assigns new values to x, y, z components corresponding to argument's components values.
        */
        pure nothrow @safe @property void xyz(Vector3 v)    { x = v.x; y = v.y; z = v.z; }        

        const string toString() { return format("[",x,", ",y,", ", z, ", ", w, "]"); }

    }
    
    /** Returns: Dot product between passed vectors. */
    pure nothrow @safe real dot(Vector4 a, Vector4 b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    }
    
    /** Returns: Outer product between passed vectors. */
    pure nothrow @safe Matrix44 outer(Vector4 a, Vector4 b)
    {
        return Matrix44( a.x * b.x,  a.x * b.y,  a.x * b.z, a.x * b.w,
                         a.y * b.x,  a.y * b.y,  a.y * b.z, a.y * b.w,
                         a.z * b.x,  a.z * b.y,  a.z * b.z, a.z * b.w,
                         a.w * b.x,  a.w * b.y,  a.w * b.z, a.w * b.w);
    }

    alias EqualityByNorm!(Vector4).equal equal; /// Introduces approximate equality function for Vector4.
    alias Lerp!(Vector4).lerp lerp;             /// Introduces linear interpolation function for Vector4.
    
    /************************************************************************************
    _Quaternion.

    For formal definition of quaternion, meaning of possible operations and related
    information see $(LINK http://en.wikipedia.org/wiki/_Quaternion),
    $(LINK http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation).
    *************************************************************************************/
    struct Quaternion
    {
        align(1)
        {
            union
            {
                struct
                {
                    float_t x; /// Components of quaternion.
                    float_t y; /// ditto
                    float_t z; /// ditto
                }
                
                Vector3 vector; /// Vector part. Can be used instead of explicit members x, y and z.
            }
    
            union
            {
                float_t w;      /// Scalar part.
                float_t scalar; /// Another name for _scalar part.
            }
        }
        
        /// Identity quaternion.
        static Quaternion identity;
        /// Quaternion with all components set to NaN.
        static Quaternion nan = { x: float_t.nan, y: float_t.nan, z: float_t.nan, w: float_t.nan };
    
        /**
        Methods to construct quaternion in C-like syntax.

        Examples:
        ------------
        Quaternion q1 = Quaternion(0, 0, 0, 1);
        Quaternion q2 = Quaternion(Vector3(0, 0, 0), 1);
        Quaternion q3 = Quaternion(Matrix33.rotationY(PI / 6), 1);
        ------------
        */
        pure nothrow @safe static Quaternion opCall(float_t x, float_t y, float_t z, float_t w)
        {
            Quaternion q;
            q.set(x, y, z, w);
            return q;
        }
    
        /** ditto */
        pure nothrow @safe static Quaternion opCall(Vector3 vector, float_t scalar)
        {
            Quaternion q;
            q.set(vector, scalar);
            return q;
        }
        
        /** ditto */
        pure @safe static Quaternion opCall(Matrix33 mat)
        {
            Quaternion q;
            q.set(mat);
            return q;
        }
        
        /** Sets values of components according to passed values. */
        pure nothrow @safe void set(float_t x, float_t y, float_t z, float_t w)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
        }
    
        /** ditto */
        pure nothrow @safe void set(Vector3 vector, float_t scalar)
        {
            this.vector = vector;
            this.scalar = scalar;
        }
        
        /**
        Sets quaternion, so that, it will represent same rotation as in mat matrix argument.
        Params:
            mat = Matrix to extract rotation from. Should be pure rotation matrix.
        Throws:
            AssertError if mat is not pure rotation matrix and module was compiled with
            contract checkings are enabled.
        */
        // NOTE: refactor to use mat.ptr instead of [] operator if
        // perforance will be unsatisfactory.
        pure @trusted void set(Matrix33 mat)
        in { assert(mat.isRotation()); }
        body
        {
            // Algorithm stolen from OGRE (http://ogre.sourceforge.net)
            real trace = mat[0, 0] + mat[1, 1] + mat[2, 2];
            real root;
        
            if ( trace > 0 )
            {
                // |w| > 1/2, may as well choose w > 1/2
                root = sqrt(trace + 1);  // 2w
                w = 0.5 * root;
                root = 0.5 / root;  // 1/(4w)
                x = (mat[2, 1] - mat[1, 2]) * root;
                y = (mat[0, 2] - mat[2, 0]) * root;
                z = (mat[1, 0] - mat[0, 1]) * root;
            }
            else
            {
                // |w| <= 1/2
                static immutable int[3] next = [ 1, 2, 0 ];
                int i = 0;
                if ( mat[1, 1] > mat[0, 0] )
                    i = 1;
                if ( mat[2, 2] > mat[i, i] )
                    i = 2;
                int j = next[i];
                int k = next[j];
                
                root = sqrt(mat[i, i] - mat[j, j] - mat[k, k] + 1);
                *(&x + i) = 0.5 * root;
                root = 0.5 / root;
                /+
                 // User "Helmut Duregger reports this sometimes 
                 // causes mirroring of rotations, and that Ogre
                 // actually uses the uncommented version below.

                 w = (mat[j, k] - mat[k, j]) * root;
                 *(&x + j) = (mat[i, j] + mat[j, i]) * root;
                 *(&x + k) = (mat[i, k] + mat[k, i]) * root;
                +/
                w = (mat[k, j] - mat[j, k]) * root;
                *(&x + j) = (mat[j, i] + mat[i, j]) * root;
                *(&x + k) = (mat[k, i] + mat[i, k]) * root;

            }
        }
        
        /** Construct quaternion that represents rotation around corresponding axis. */
        pure nothrow @safe static Quaternion rotationX(float_t radians)
        {
            Quaternion q;
            
            float_t s = sin(radians * 0.5f);
            float_t c = cos(radians * 0.5f);
            q.x = s;
            q.y = 0;
            q.z = 0;
            q.w = c;
            
            return q;
        }
    
        /** ditto */
        pure nothrow @safe static Quaternion rotationY(float_t radians)
        {
            Quaternion q;
            
            float_t s = sin(radians * 0.5f);
            float_t c = cos(radians * 0.5f);
            q.x = 0;
            q.y = s;
            q.z = 0;
            q.w = c;
            
            return q;
        }
    
        /** ditto */
        pure nothrow @safe static Quaternion rotationZ(float_t radians)
        {
            Quaternion q;
            
            float_t s = sin(radians * 0.5f);
            float_t c = cos(radians * 0.5f);
            q.x = 0;
            q.y = 0;
            q.z = s;
            q.w = c;
            
            return q;
        }
    
        /**
        Constructs quaternion that represents _rotation specified by euler angles passed as arguments.
        Order of _rotation application is: roll (Z axis), pitch (X axis), yaw (Y axis).
        */
        pure nothrow @safe static Quaternion rotation(float_t yaw, float_t pitch, float_t roll)
        {
            return Quaternion.rotationY(yaw) * Quaternion.rotationX(pitch) * Quaternion.rotationZ(roll);
        }
    
        /**
        Constructs quaternion that represents _rotation around 'axis' _axis by 'radians' angle.
        */
        pure nothrow @safe static Quaternion rotation(Vector3 axis, float_t radians)
        {
            Quaternion q;
            
            float_t s = sin(radians * 0.5f);
            float_t c = cos(radians * 0.5f);
            q.x = axis.x * s;
            q.y = axis.y * s;
            q.z = axis.z * s;
            q.w = c;
            
            return q;
        }
    
        /** Returns: Norm (also known as length, magnitude) of quaternion. */
        const pure nothrow @safe @property real norm()
        {
            return sqrt(x*x + y*y + z*z + w*w);
        }
    
        /**
        Returns: Square of vector's norm.
        
        Method doesn't need calculation of square root.
        */
        const pure nothrow @safe @property real normSquare()
        {
            return x*x + y*y + z*z + w*w;
        }
    
        /** Normalizes this quaternion. Returns: the original length. */
        pure @safe real normalize()
        {
            float_t n = norm();
            assert( greater(n, 0) );
            this /= n;
            return n;
        }
    
        /** Returns: Normalized copy of this quaternion. */
        const pure @safe Quaternion normalized()
        {
            float_t n = norm();
            assert( greater(n, 0) );
            return Quaternion(x / n, y / n, z / n, w / n);
        }
    
        /**
        Returns: Whether this quaternion is unit.
        Params:
            relprec, absprec = Parameters passed to equal function while comparison of
                               norm square and 1. Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isUnit(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal( normSquare(), 1, relprec, absprec );
        }
    
        /** Returns: Conjugate quaternion. */
        const pure nothrow @safe @property Quaternion conj()
        {
            return Quaternion(-vector, scalar);
        }
    
        /** Invert this quaternion. */
        pure @safe void invert()
        {
            float_t n = norm();
            assert( greater(n, 0) );
            vector = -vector / n;
            scalar =  scalar / n;
        }
    
        /** Returns: Inverse copy of this quaternion. */
        const pure @safe @property Quaternion inverse()
        {
            float_t n = norm();
            assert( greater(n, 0) );
            return conj / n;
        }
        
        /**
        Returns: Extracted euler angle with the assumption that rotation is applied in order:
        _roll (Z axis), _pitch (X axis), _yaw (Y axis).
        */
        const pure nothrow @safe @property real yaw()
        {
            return atan2(2 * (w*y + x*z), w*w - x*x - y*y + z*z);
        }
    
        /** ditto */
        const pure nothrow @safe @property real pitch()
        {
            return asin(2 * (w*x - y*z));
        }
        
        /** ditto */
        const pure nothrow @safe @property real roll()
        {
            return atan2(2 * (w*z + x*y), w*w - x*x + y*y - z*z);
        }
        
        /** Returns: Whether all components are normalized. */
        const pure nothrow @safe bool isNormal()
        {
            return std.math.isNormal(x) && std.math.isNormal(y) && std.math.isNormal(z) && std.math.isNormal(w);
        }
    
        /** Returns: float_t pointer to x component of this vector. It's like a _ptr method for arrays. */
        pure nothrow @safe float_t* ptr()
        {
            return cast(float_t*)(&x);
        }
    
        /** Returns: Component corresponded to passed index. */
        pure @safe float_t opIndex(size_t ort)
        in { assert(ort <= Ort.W); }
        body
        {
            return (cast(float_t*)&this)[cast(int)ort];
        }
    
        /** Assigns new _value to component corresponded to passed index. */
        pure @safe void opIndexAssign(float_t value, size_t ort)
        in { assert(ort <= Ort.W); }
        body
        {
            (cast(float_t*)&this)[cast(int)ort] = value;
        }
    
        /**
        Standard operators that have meaning exactly the same as for Vector4.
        
        Note that division operators do no cheks of value of k, so in case of division
        by 0 result vector will have infinity components. You can check this with isNormal()
        method.
        */
        const pure nothrow @safe bool opEquals(Quaternion q)
        {
            return x == q.x && y == q.y && z == q.z && w == q.w;
        }
    
        /** ditto */
        const pure nothrow @safe Quaternion opNeg()
        {
            return Quaternion(-x, -y, -z, -w);
        }
    
        /** ditto */
        const pure nothrow @safe Quaternion opAdd(Quaternion q)
        {
            return Quaternion(x + q.x, y + q.y, z + q.z, w + q.w);
        }
    
        /** ditto */
        pure nothrow @safe void opAddAssign(Quaternion q)
        {
            x += q.x;
            y += q.y;
            z += q.z;
            w += q.w;
        }
    
        /** ditto */
        const pure nothrow @safe Quaternion opSub(Quaternion q)
        {
            return Quaternion(x - q.x, y - q.y, z - q.z, w - q.w);
        }
    
        /** ditto */
        pure nothrow @safe void opSubAssign(Quaternion q)
        {
            x -= q.x;
            y -= q.y;
            z -= q.z;
            w -= q.w;
        }
    
        /** ditto */
        const pure nothrow @safe Quaternion opMul(float_t k)
        {
            return Quaternion(x * k, y * k, z * k, w * k);
        }

        /** ditto */
        const pure nothrow @safe Quaternion opMul_r(float_t k)
        {
            return Quaternion(x * k, y * k, z * k, w * k);
        }
    
        /** ditto */
        const pure nothrow @safe Quaternion opDiv(float_t k)
        {
            return Quaternion(x / k, y / k, z / k, w / k);
        }
    
        /** ditto */
        pure nothrow @safe void opDivAssign(float_t k)
        {
            x /= k;
            y /= k;
            z /= k;
            w /= k;
        }
    
        /**
        Quaternion multiplication operators. Result of Q1*Q2 is quaternion that represents
        rotation which has meaning of application Q2's rotation and the Q1's rotation.
        */
        const pure nothrow @safe Quaternion opMul(Quaternion q)
        {
            return Quaternion(
                w * q.x + x * q.w + y * q.z - z * q.y,
                w * q.y + y * q.w + z * q.x - x * q.z,
                w * q.z + z * q.w + x * q.y - y * q.x,
                w * q.w - x * q.x - y * q.y - z * q.z );
        }
        
        /** ditto */
        pure nothrow @safe void opMulAssign(Quaternion q)
        {
            set(w * q.x + x * q.w + y * q.z - z * q.y,
                w * q.y + y * q.w + z * q.x - x * q.z,
                w * q.z + z * q.w + x * q.y - y * q.x,
                w * q.w - x * q.x - y * q.y - z * q.z );
        }
    
        /** Returns: Copy of this quaternion with float type components. */
        const pure nothrow @safe Quaternionf toQuaternionf()
        {
            return Quaternionf(cast(float)x, cast(float)y, cast(float)z, cast(float)w);
        }
        
        /** Returns: Copy of this vector with double type components. */
        const pure nothrow @safe Quaterniond toQuaterniond()
        {
            return Quaterniond(cast(double)x, cast(double)y, cast(double)z, cast(double)w);
        }
        
        /** Returns: Copy of this vector with real type components. */
        const pure nothrow @safe Quaternionr toQuaternionr()
        {
            return Quaternionr(cast(real)x, cast(real)y, cast(real)z, cast(real)w);
        }

        const string toString() { return format("[",x,", ",y,", ", z, ", ", w, "]"); }
    }    
    
    alias EqualityByNorm!(Quaternion).equal equal;  /// Introduces approximate equality function for Quaternion.
    alias Lerp!(Quaternion).lerp lerp;              /// Introduces linear interpolation function for Quaternion.
    
    /**
    Returns:
        Product of spherical linear interpolation between q0 and q1 with parameter t.
    References:
        $(LINK http://en.wikipedia.org/wiki/Slerp).
    */
    @trusted Quaternion slerp(Quaternion q0, Quaternion q1, real t)
    {
        real cosTheta = q0.x * q1.x + q0.y * q1.y + q0.z * q1.z + q0.w * q1.w;
        real theta = acos(cosTheta);
    
        if ( equal(fabs(theta), 0) )
            return lerp(q0, q1, t);
    
        real sinTheta = sin(theta);
        real coeff0 = sin((1 - t) * theta) / sinTheta;
        real coeff1 = sin(t * theta) / sinTheta;
        
        // Invert rotation if necessary
        if (cosTheta < 0.0f)
        {
            coeff0 = -coeff0;
            // taking the complement requires renormalisation
            Quaternion ret = coeff0 * q0 + coeff1 * q1;
            return ret.normalized();
        }
        
        return coeff0 * q0 + coeff1 * q1;    
    }
    
    /************************************************************************************
    2x2 Matrix.

    $(LINK http://en.wikipedia.org/wiki/Transformation_matrix).
    *************************************************************************************/
    struct Matrix22
    {
        align(1) union
        {
            struct
            {
                float_t m00, m10;
                float_t m01, m11;
            }
    
            float_t[2][2] m;
            Vector2[2]    v;
            float_t[4]    a;
        }
    
        /// Identity matrix.
        static immutable Matrix22 identity = {
            1, 0,
            0, 1 };
        /// Matrix with all elements set to NaN.
        static immutable Matrix22 nan = {
            float_t.nan, float_t.nan,
            float_t.nan, float_t.nan, };
        /// Matrix with all elements set to 0.
        static immutable Matrix22 zero = {
            0, 0,
            0, 0 };
    
        /**
        Methods to construct matrix in C-like syntax.

        In case with array remember about column-major matrix memory layout,
        note last line with assert in example.

        Examples:
        ------------
        Matrix22 mat1 = Matrix22(1,2,3,4);
        static float[9] a = [1,2,3,4];
        Matrix22 mat2 = Matrix22(a);

        assert(mat1 == mat2.transposed);
        ------------
        */
        pure nothrow @safe static Matrix22 opCall(float_t m00, float_t m01,
                               float_t m10, float_t m11)
        {
            Matrix22 mat;
            mat.m00 = m00;        mat.m01 = m01;
            mat.m10 = m10;        mat.m11 = m11;
            return mat;
        }
        
        /** ditto */
        pure @safe static Matrix22 opCall(float_t[4] a)
        {
            Matrix22 mat;
            mat.a[0..4] = a[0..4].dup;
            return mat;
        }
        
        /**
        Method to construct matrix in C-like syntax. Sets columns to passed vector
        arguments.
        */
        pure nothrow @safe static Matrix22 opCall(Vector2 basisX, Vector2 basisY)
        {
            Matrix22 mat;
            mat.v[0] = basisX;
            mat.v[1] = basisY;
            return mat;
        }
    
        /** Sets elements to passed values. */
        pure nothrow @safe void set(float_t m00, float_t m01,
                 float_t m10, float_t m11)
        {
            this.m00 = m00;        this.m01 = m01;
            this.m10 = m10;        this.m11 = m11;
        }
        
        /** Sets elements as _a copy of a contents. Remember about column-major matrix memory layout. */
        pure @safe void set(float_t[4] a)
        {
            this.a[0..4] = a[0..4].dup;
        }
        
        /** Sets columns to passed basis vectors. */
        pure nothrow @safe void set(Vector2 basisX, Vector2 basisY)
        {
            v[0] = basisX;
            v[1] = basisY;
        }
        
        /** Returns: Whether all components are normalized numbers. */
        const pure nothrow @safe bool isNormal()
        {
            return
                std.math.isNormal(m00) && std.math.isNormal(m01) &&
                std.math.isNormal(m10) && std.math.isNormal(m11);
        }
        
        /**
        Returns: Whether this matrix is identity.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                               Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isIdentity(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal(this, identity, relprec, absprec);
        }
        
        /**
        Returns: Whether this matrix is zero.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                               Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isZero(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal(normSquare, 0, relprec, absprec);
        }
        
        /**
        Returns: Whether this matrix is orthogonal.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                               Have the same meaning as in equal function.
        References:
            $(LINK http://en.wikipedia.org/wiki/Orthogonal_matrix).
        */
        const pure nothrow @safe bool isOrthogonal(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal(abs(cross(v[0],v[1])), 1.0, relprec, absprec);
        }
        
        /**
        Returns: Whether this matrix represents pure rotation. I.e. hasn't scale admixture.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                               Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isRotation(int relprec = defrelprec, int absprec = defabsprec)
        {
            return isOrthogonal(relprec, absprec);
        }
    
        /** Constructs _scale matrix with _scale coefficients specified as arguments. */
        pure nothrow @safe static Matrix22 scale(float_t x, float_t y)
        {
            Matrix22 mat = identity;
            with (mat)
            {
                m00 = x;
                m11 = y;
            }
    
            return mat;
        }
    
        /** Constructs _scale matrix with _scale coefficients specified as v's components. */
        pure nothrow @safe static Matrix22 scale(Vector2 v)
        {
            return scale(v.x, v.y);
        }
    
        /** Construct matrix that represents rotation. */
        pure nothrow @safe static Matrix22 rotation(float_t radians)
        {
            Matrix22 mat = identity;
            float_t c = cos(radians);
            float_t s = sin(radians);
            with (mat)
            {
                m00 = m11 = c;
                m10 = s;
                m01 = -s;
            }
    
            return mat;
        }
    
    
        /**
        Constructs matrix that represents _rotation same as in passed in complex number q.
        Method works with assumption that q is unit.
        Throws:
            AssertError on non-unit quaternion call attempt if you compile with
            contract checks enabled.
        */
/*
        static Matrix22 rotation(complex q)
        in { assert( q.isUnit() ); }
        body
        {
            float_t tx  = 2.f * q.x;
            float_t ty  = 2.f * q.y;
            float_t tz  = 2.f * q.z;
            float_t twx = tx * q.w;
            float_t twy = ty * q.w;
            float_t twz = tz * q.w;
            float_t txx = tx * q.x;
            float_t txy = ty * q.x;
            float_t txz = tz * q.x;
            float_t tyy = ty * q.y;
            float_t tyz = tz * q.y;
            float_t tzz = tz * q.z;
            
            Matrix22 mat;
            with (mat)
            {
                m00 = 1.f - (tyy + tzz);    m01 = txy + twz;            m02 = txz - twy;        
                m10 = txy - twz;            m11 = 1.f - (txx + tzz);    m12 = tyz + twx;        
                m20 = txz + twy;            m21 = tyz - twx;            m22 = 1.f - (txx + tyy);
            }
            
            return mat;
        }
*/        
        /**
        Returns: Inverse copy of this matrix.
        
        In case if this matrix is singular (i.e. determinant = 0) result matrix will has
        infinity elements. You can check this with isNormal() method.
        */
        const pure nothrow @safe @property Matrix22 inverse()
        {
            Matrix22 mat;
            
            mat.m00 =  m11;
            mat.m01 = -m01;
            mat.m10 = -m10;
            mat.m11 =  m00;

            real det = m00 * m11 - m01 * m10;
            
            for (int i = 4; i--; )
                mat.a[i] /= det;
    
            return mat;
        }
        
        /**
        Inverts this matrix.
        
        In case if matrix is singular (i.e. determinant = 0) result matrix will has
        infinity elements. You can check this with isNormal() method.
        */
        pure nothrow @safe void invert()
        {
            real idet = 1.0/(m00 * m11 - m01 * m10);
            swap(m00,m11);
            m10 = -m10;
            m01 = -m01;
            this *= idet;
        }
        
        /** Returns: Determinant */
        const pure nothrow @safe @property real determinant()
        {
            return m00 * m11 - m10 * m01;
        }
        
        /**
        Returns: Frobenius _norm of matrix. I.e. square root from sum of all elements' squares.
        References:
            $(LINK http://en.wikipedia.org/wiki/Frobenius_norm#Frobenius_norm).
        */
        const pure nothrow @safe @property real norm()
        {
            return sqrt( normSquare );
        }
        
        /**
        Returns: Square of Frobenius _norm of matrix. I.e. sum of all elements' squares.

        Method doesn't need calculation of square root.

        References:
            $(LINK http://en.wikipedia.org/wiki/Frobenius_norm#Frobenius_norm).
        */
        const pure nothrow @safe @property real normSquare()
        {
            real ret = 0;
            for (int i = 4; i--; )
            {
                real x = a[i];
                ret += x * x;
            }
            
            return ret;
        }
        
        /** Transposes this matrix. */
        pure nothrow @safe void transpose()
        {
            /*           */        swap(m01, m10);
            /*           */        /*           */
        }
        
        /** Returns: Transposed copy of this matrix. */
        const pure nothrow @safe @property Matrix22 transposed()
        {
            return Matrix22(
                m00, m10,
                m01, m11 );
        }
        
        /**
        Makes polar decomposition of this matrix. Denote this matrix with 'M', in
        that case M=Q*S.
        
        Method is useful to decompose your matrix into rotation 'Q' and scale+shear 'S'. If you
        didn't use shear transform matrix S will be diagonal, i.e. represent scale. This can
        have many applications, particulary you can use method for suppressing errors in pure
        rotation matrices after long multiplication chain.
        
        Params:
            Q = Output matrix, will be orthogonal after decomposition.
                Argument shouldn't be null.
            S = Output matrix, will be symmetric non-negative definite after
                decomposition. Argument shouldn't be null.

        Examples:
        --------
        Matrix22 Q, S;
        Matrix22 rot = Matrix22.rotationZ(PI / 7);
        Matrix22 scale = Matrix22.scale(-1, 2, 3);
        Matrix22 composition = rot * scale;
        composition.polarDecomposition(Q, S);    
        assert( equal(Q * S, composition) );
        --------

        References:
            $(LINK http://www.cs.wisc.edu/graphics/Courses/cs-838-2002/Papers/polar-decomp.pdf)
        */
        const @trusted void polarDecomposition(out Matrix22 Q, out Matrix22 S)
            out { assert(Q.isRotation(), 
                         "(postcondition) Q not a rotation:\n" ~ Q.toString()); }
        body
        {
            // TODO: Optimize, we need only sign of determinant, not value
            if (determinant < 0)
                Q = this * (-identity);
            else
                Q = this;
                
            // use scaled Newton method to orthonormalize Q
            int maxIterations = 100; // avoid deadlock
            Matrix22 Qp = Q;
            Q = 0.5f * (Q + Q.transposed.inverse);
            while (!(Q - Qp).isZero() && maxIterations--)
            {
                Matrix22 Qinv = Q.inverse;
                real gamma = sqrt( Qinv.norm / Q.norm );
                Qp = Q;
                Q = 0.5f * (gamma * Q + (1 / gamma) * Qinv.transposed);
            }
            
            assert( maxIterations != -1 );
            
            S = Q.transposed * this;
        }
    
        /**
        Standard operators that have intuitive meaning, same as in classical math.
        
        Note that division operators do no cheks of value of k, so in case of division
        by 0 result matrix will have infinity components. You can check this with isNormal()
        method.
        */
        const pure nothrow @safe Matrix22 opNeg()
        {
            return Matrix22(-m00, -m01,
                            -m10, -m11);
        }
    
        /** ditto */
        const pure nothrow @safe Matrix22 opAdd(Matrix22 mat)
        {
            return Matrix22(m00 + mat.m00, m01 + mat.m01,
                            m10 + mat.m10, m11 + mat.m11);
        }
    
        /** ditto */
        pure nothrow @safe void opAddAssign(Matrix22 mat)
        {
            m00 += mat.m00; m01 += mat.m01;
            m10 += mat.m10; m11 += mat.m11;
        }
    
        /** ditto */
        const pure nothrow @safe Matrix22 opSub(Matrix22 mat)
        {
            return Matrix22(m00 - mat.m00, m01 - mat.m01,
                            m10 - mat.m10, m11 - mat.m11);
        }
    
        /** ditto */
        pure nothrow @safe void opSubAssign(Matrix22 mat)
        {
            m00 -= mat.m00; m01 -= mat.m01;
            m10 -= mat.m10; m11 -= mat.m11;
        }
    
        /** ditto */
        const pure nothrow @safe Matrix22 opMul(float_t k)
        {
            return Matrix22(m00 * k, m01 * k,
                            m10 * k, m11 * k);
        }
    
        /** ditto */
        pure nothrow @safe void opMulAssign(float_t k)
        {
            m00 *= k; m01 *= k;
            m10 *= k; m11 *= k;
        }
    
        /** ditto */
        const pure nothrow @safe Matrix22 opMul_r(float_t k)
        {
            return Matrix22(m00 * k, m01 * k,
                            m10 * k, m11 * k);
        }
    
        /** ditto */
        const pure nothrow @safe Matrix22 opDiv(float_t k)
        {
            return Matrix22(m00 / k, m01 / k,
                            m10 / k, m11 / k);
        }
    
        /** ditto */
        pure nothrow @safe void opDivAssign(float_t k)
        {
            m00 /= k; m01 /= k;
            m10 /= k; m11 /= k;
        }
    
        /** ditto */
        const pure nothrow @safe bool opEquals(Matrix22 mat)
        {
            return m00 == mat.m00 && m01 == mat.m01 &&
                   m10 == mat.m10 && m11 == mat.m11;
        }

        /** ditto */
        const pure nothrow @safe Matrix22 opMul(Matrix22 mat)
        {
            return Matrix22(
                m00 * mat.m00 + m01 * mat.m10,   m00 * mat.m01 + m01 * mat.m11,
                m10 * mat.m00 + m11 * mat.m10,   m10 * mat.m01 + m11 * mat.m11 );
        }
    
        /** ditto */
        pure nothrow @safe void opMulAssign(Matrix22 mat)
        {
            this = this * mat;
        }
    
        /** ditto */
        const pure nothrow @safe Vector2 opMul(Vector2 v)
        {
            return Vector2(v.x * m00 + v.y * m01,
                           v.x * m10 + v.y * m11 );
        }
    
        /** Returns: Element at row'th _row and col'th column. */
        pure @safe float_t opIndex(uint row, uint col)
        in { assert( row < 2 && col < 2 ); }
        body
        {
            return m[col][row];
        }
    
        /** Assigns value f to element at row'th _row and col'th column. */
        pure @safe void opIndexAssign(float_t f, uint row, uint col)
        in { assert( row < 2 && col < 2 ); }
        body
        {
            m[col][row] = f;
        }
        
        /** Returns: Vector representing col'th column. */
        pure @safe Vector2 opIndex(uint col)
        in { assert( col < 2 ); }
        body
        {
            return v[col];
        }
        
        /** Replaces elements in col'th column with v's values. */
        pure @safe Vector2 opIndexAssign(Vector2 v, uint col)
        in { assert( col < 2 ); }
        body
        {
            return this.v[col] = v;
        }
    
        /**
        Returns: float_t pointer to [0,0] element of this matrix. It's like a _ptr method for arrays.
        
        Remember about column-major matrix memory layout.
        */
        pure nothrow @safe float_t* ptr()
        {
            return a.ptr;
        }
        
        /** Returns: Copy of this matrix with float type elements. */
        const pure nothrow @safe Matrix22f toMatrix22f()
        {
            return Matrix22f(
                cast(float)m00, cast(float)m01,
                cast(float)m10, cast(float)m11 );
        }
        
        /** Returns: Copy of this matrix with double type elements. */
        const pure nothrow @safe Matrix22d toMatrix22d()
        {
            return Matrix22d(
                cast(double)m00, cast(double)m01,
                cast(double)m10, cast(double)m11 );
        }
        
        /** Returns: Copy of this matrix with real type elements. */
        const pure nothrow @safe Matrix22r toMatrix22r()
        {
            return Matrix22r(
                cast(real)m00, cast(real)m01,
                cast(real)m10, cast(real)m11     );
        }

        const string toString() { 
            return format("[" ,m00, ", " ,m01, ",\n",
                          " " ,m10, ", " ,m11, "]");
        }
    }
    
    
    alias EqualityByNorm!(Matrix22).equal equal; /// Introduces approximate equality function for Matrix22.
    alias Lerp!(Matrix22).lerp lerp;             /// Introduces linear interpolation function for Matrix22.

    /************************************************************************************
    3x3 Matrix.

    For formal definition of quaternion, meaning of possible operations and related
    information see $(LINK http://en.wikipedia.org/wiki/Matrix(mathematics)),
    $(LINK http://en.wikipedia.org/wiki/Transformation_matrix).
    *************************************************************************************/
    struct Matrix33
    {
        align(1) union
        {
            struct
            {
                float_t m00, m10, m20;
                float_t m01, m11, m21;
                float_t m02, m12, m22;
            }
    
            float_t[3][3] m;
            Vector3[3]    v;
            float_t[9]    a;
        }
    
        /// Identity matrix.
        static immutable Matrix33 identity = {
            1, 0, 0,
            0, 1, 0,
            0, 0, 1 };
        /// Matrix with all elements set to NaN.
        static immutable Matrix33 nan = {
            float_t.nan, float_t.nan, float_t.nan,
            float_t.nan, float_t.nan, float_t.nan,
            float_t.nan, float_t.nan, float_t.nan };
        /// Matrix with all elements set to 0.
        static immutable Matrix33 zero = {
            0, 0, 0,
            0, 0, 0,
            0, 0, 0 };
    
        /**
        Methods to construct matrix in C-like syntax.

        In case with array remember about column-major matrix memory layout,
        note last line with assert in example.

        Examples:
        ------------
        Matrix33 mat1 = Matrix33(1,2,3,4,5,6,7,8,9);
        static float[9] a = [1,2,3,4,5,6,7,8,9];
        Matrix33 mat2 = Matrix33(a);

        assert(mat1 == mat2.transposed);
        ------------
        */
        pure nothrow @safe static Matrix33 opCall(float_t m00, float_t m01, float_t m02,
                               float_t m10, float_t m11, float_t m12,
                               float_t m20, float_t m21, float_t m22)
        {
            Matrix33 mat;
            mat.m00 = m00;        mat.m01 = m01;        mat.m02 = m02;
            mat.m10 = m10;        mat.m11 = m11;        mat.m12 = m12;
            mat.m20 = m20;        mat.m21 = m21;        mat.m22 = m22;
            return mat;
        }
        
        /** ditto */
        pure @safe static Matrix33 opCall(float_t[9] a)
        {
            Matrix33 mat;
            mat.a[0..9] = a[0..9].dup;
            return mat;
        }
        
        /**
        Method to construct matrix in C-like syntax. Sets columns to passed vector
        arguments.
        */
        pure nothrow @safe static Matrix33 opCall(Vector3 basisX, Vector3 basisY, Vector3 basisZ)
        {
            Matrix33 mat;
            mat.v[0] = basisX;
            mat.v[1] = basisY;
            mat.v[2] = basisZ;
            return mat;
        }
    
        /** Sets elements to passed values. */
        pure nothrow @safe void set(float_t m00, float_t m01, float_t m02,
                 float_t m10, float_t m11, float_t m12,
                 float_t m20, float_t m21, float_t m22)
        {
            this.m00 = m00;        this.m01 = m01;        this.m02 = m02;
            this.m10 = m10;        this.m11 = m11;        this.m12 = m12;
            this.m20 = m20;        this.m21 = m21;        this.m22 = m22;
        }
        
        /** Sets elements as _a copy of a contents. Remember about column-major matrix memory layout. */
        pure @safe void set(float_t[9] a)
        {
            this.a[0..9] = a[0..9].dup;
        }
        
        /** Sets columns to passed basis vectors. */
        pure nothrow @safe void set(Vector3 basisX, Vector3 basisY, Vector3 basisZ)
        {
            v[0] = basisX;
            v[1] = basisY;
            v[2] = basisZ;
        }
        
        /** Returns: Whether all components are normalized numbers. */
        const pure nothrow @safe bool isNormal()
        {
            return
                std.math.isNormal(m00) && std.math.isNormal(m01) && std.math.isNormal(m02) &&
                std.math.isNormal(m10) && std.math.isNormal(m11) && std.math.isNormal(m12) &&
                std.math.isNormal(m20) && std.math.isNormal(m21) && std.math.isNormal(m22);
        }
        
        /**
        Returns: Whether this matrix is identity.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                               Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isIdentity(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal(this, identity, relprec, absprec);
        }
        
        /**
        Returns: Whether this matrix is zero.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                               Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isZero(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal(normSquare(), 0, relprec, absprec);
        }
        
        /**
        Returns: Whether this matrix is orthogonal.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                               Have the same meaning as in equal function.
        References:
            $(LINK http://en.wikipedia.org/wiki/Orthogonal_matrix).
        */
        const pure nothrow @safe bool isOrthogonal(int relprec = defrelprec, int absprec = defabsprec)
        {
            return isBasisOrthonormal(v[0], v[1], v[2], relprec, absprec);
        }
        
        /**
        Returns: Whether this matrix represents pure rotation. I.e. hasn't scale admixture.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                               Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isRotation(int relprec = defrelprec, int absprec = defabsprec)
        {
            return isOrthogonal(relprec, absprec) && equal(v[2], cross(v[0], v[1]), relprec, absprec);
        }
    
        /** Constructs _scale matrix with _scale coefficients specified as arguments. */
        pure nothrow @safe static Matrix33 scale(float_t x, float_t y, float_t z)
        {
            Matrix33 mat = identity;
            with (mat)
            {
                m00 = x;
                m11 = y;
                m22 = z;
            }
    
            return mat;
        }
    
        /** Constructs _scale matrix with _scale coefficients specified as v's components. */
        pure nothrow @safe static Matrix33 scale(Vector3 v)
        {
            return scale(v.x, v.y, v.z);
        }
    
        /** Construct matrix that represents rotation around corresponding axis. */
        pure nothrow @safe static Matrix33 rotationX(float_t radians)
        {
            Matrix33 mat = identity;
            float_t c = cos(radians);
            float_t s = sin(radians);
            with (mat)
            {
                m11 = m22 = c;
                m21 = s;
                m12 = -s;            
            }
    
            return mat;
        }
    
        /** ditto */
        pure nothrow @safe static Matrix33 rotationY(float_t radians)
        {
            Matrix33 mat = identity;
            float_t c = cos(radians);
            float_t s = sin(radians);
            with (mat)
            {
                m00 = m22 = c;
                m20 = -s;
                m02 = s;            
            }
    
            return mat;
        }
    
        /** ditto */
        pure nothrow @safe static Matrix33 rotationZ(float_t radians)
        {
            Matrix33 mat = identity;
            float_t c = cos(radians);
            float_t s = sin(radians);
            with (mat)
            {
                m00 = m11 = c;
                m10 = s;
                m01 = -s;            
            }
    
            return mat;
        }
    
        /**
        Constructs matrix that represents _rotation specified by euler angles passed as arguments.
        Order of _rotation application is: roll (Z axis), pitch (X axis), yaw (Y axis).
        */
        pure nothrow @safe static Matrix33 rotation(float_t yaw, float_t pitch, float_t roll)
        {
            return Matrix33.rotationY(yaw) * Matrix33.rotationX(pitch) * Matrix33.rotationZ(roll);
        }
    
        /**
        Constructs matrix that represents _rotation specified by axis and angle.
        Method works with assumption that axis is unit vector.        
        Throws:
            AssertError on non-unit axis call attempt if module was compiled with
            contract checks enabled.
        */
        pure @safe static Matrix33 rotation(Vector3 axis, float_t radians)
        in { assert( axis.isUnit() ); }
        body
        {
            real c = cos(radians);
            real s = sin(radians);
            real cc = 1.0 - c;
            real x2 = axis.x * axis.x;
            real y2 = axis.y * axis.y;
            real z2 = axis.z * axis.z;
            real xycc = axis.x * axis.y * cc;
            real xzcc = axis.x * axis.z * cc;
            real yzcc = axis.y * axis.z * cc;
            real xs = axis.x * s;
            real ys = axis.y * s;
            real zs = axis.z * s;
    
            Matrix33 mat;
            with (mat)
            {
                m00 = x2 * cc + c;      m01 = xycc - zs;        m02 = xzcc + ys;
                m10 = xycc + zs;        m11 = y2 * cc + c;      m12 = yzcc - xs;
                m20 = xzcc - ys;        m21 = yzcc + xs;        m22 = z2 * cc + c;
            }
    
            return mat;
        }
        
        /**
        Constructs matrix that represents _rotation same as in passed quaternion q.
        Method works with assumption that q is unit.
        Throws:
            AssertError on non-unit quaternion call attempt if you compile with
            contract checks enabled.
        */
        pure @safe static Matrix33 rotation(Quaternion q)
        in { assert( q.isUnit() ); }
        body
        {
            float_t tx  = 2.0f * q.x;
            float_t ty  = 2.0f * q.y;
            float_t tz  = 2.0f * q.z;
            float_t twx = tx * q.w;
            float_t twy = ty * q.w;
            float_t twz = tz * q.w;
            float_t txx = tx * q.x;
            float_t txy = ty * q.x;
            float_t txz = tz * q.x;
            float_t tyy = ty * q.y;
            float_t tyz = tz * q.y;
            float_t tzz = tz * q.z;
            
            Matrix33 mat;
            with (mat)
            {
                m00 = 1.0f - (tyy + tzz);   m01 = txy - twz;            m02 = txz + twy;
                m10 = txy + twz;            m11 = 1.0f - (txx + tzz);   m12 = tyz - twx;
                m20 = txz - twy;            m21 = tyz + twx;            m22 = 1.0f - (txx + tyy);
            }
            
            return mat;
        }
        
        /**
        Returns: Inverse copy of this matrix.
        
        In case if this matrix is singular (i.e. determinant = 0) result matrix will has
        infinity elements. You can check this with isNormal() method.
        */
        const pure nothrow @safe @property Matrix33 inverse()
        {
            Matrix33 mat;
            
            mat.m00 = m11 * m22 - m12 * m21;
            mat.m01 = m02 * m21 - m01 * m22;
            mat.m02 = m01 * m12 - m02 * m11;
            mat.m10 = m12 * m20 - m10 * m22;
            mat.m11 = m00 * m22 - m02 * m20;
            mat.m12 = m02 * m10 - m00 * m12;
            mat.m20 = m10 * m21 - m11 * m20;
            mat.m21 = m01 * m20 - m00 * m21;
            mat.m22 = m00 * m11 - m01 * m10;
            
            real det = m00 * mat.m00 + m01 * mat.m10 + m02 * mat.m20;
            
            for (int i = 9; i--; )
                mat.a[i] /= det;
    
            return mat;
        }
        
        /**
        Inverts this matrix.
        
        In case if matrix is singular (i.e. determinant = 0) result matrix will has
        infinity elements. You can check this with isNormal() method.
        */
        pure nothrow @safe void invert()
        {
            this = inverse();
        }
        
        /** Returns: Determinant */
        const pure nothrow @safe @property real determinant()
        {
            real cofactor00 = m11 * m22 - m12 * m21;
            real cofactor10 = m12 * m20 - m10 * m22;
            real cofactor20 = m10 * m21 - m11 * m20;
            
            return m00 * cofactor00 + m01 * cofactor10 + m02 * cofactor20;
        }
        
        /**
        Returns: Frobenius _norm of matrix. I.e. square root from sum of all elements' squares.
        References:
            $(LINK http://en.wikipedia.org/wiki/Frobenius_norm#Frobenius_norm).
        */
        const pure nothrow @safe @property real norm()
        {
            return sqrt( normSquare );
        }
        
        /**
        Returns: Square of Frobenius _norm of matrix. I.e. sum of all elements' squares.

        Method doesn't need calculation of square root.

        References:
            $(LINK http://en.wikipedia.org/wiki/Frobenius_norm#Frobenius_norm).
        */
        const pure nothrow @safe @property real normSquare()
        {
            real ret = 0;
            for (int i = 9; i--; )
            {
                real x = a[i];
                ret += x * x;
            }
            
            return ret;
        }
        
        /** Transposes this matrix. */
        pure nothrow @safe void transpose()
        {
            /*           */        swap(m01, m10);        swap(m02, m20);
            /*           */        /*           */        swap(m12, m21);
            /*           */        /*           */        /*           */
        }
        
        /** Returns: Transposed copy of this matrix. */
        const pure nothrow @safe @property Matrix33 transposed()
        {
            return Matrix33(
                m00, m10, m20,
                m01, m11, m21,
                m02, m12, m22 );
        }
        
        /**
        Makes polar decomposition of this matrix. Denote this matrix with 'M', in
        that case M=Q*S.
        
        Method is useful to decompose your matrix into rotation 'Q' and scale+shear 'S'. If you
        didn't use shear transform matrix S will be diagonal, i.e. represent scale. This can
        have many applications, particulary you can use method for suppressing errors in pure
        rotation matrices after long multiplication chain.
        
        Params:
            Q = Output matrix, will be orthogonal after decomposition.
                Argument shouldn't be null.
            S = Output matrix, will be symmetric non-negative definite after
                decomposition. Argument shouldn't be null.

        Examples:
        --------
        Matrix33 Q, S;
        Matrix33 rot = Matrix33.rotationZ(PI / 7);
        Matrix33 scale = Matrix33.scale(-1, 2, 3);
        Matrix33 composition = rot * scale;
        composition.polarDecomposition(Q, S);    
        assert( equal(Q * S, composition) );
        --------

        References:
            $(LINK http://www.cs.wisc.edu/graphics/Courses/cs-838-2002/Papers/polar-decomp.pdf)
        */
        const pure @safe void polarDecomposition(out Matrix33 Q, out Matrix33 S)
        out { assert(Q.isRotation()); }
        body
        {
            // TODO: Optimize, we need only sign of determinant, not value
            if (determinant < 0)
                Q = this * (-identity);
            else
                Q = this;
                
            // use scaled Newton method to orthonormalize Q
            int maxIterations = 100; // avoid deadlock
            Matrix33 Qp = Q;
            Q = 0.5f * (Q + Q.transposed.inverse);
            while (!(Q - Qp).isZero() && maxIterations--)
            {
                Matrix33 Qinv = Q.inverse;
                real gamma = sqrt( Qinv.norm / Q.norm );
                Qp = Q;
                Q = 0.5f * (gamma * Q + (1 / gamma) * Qinv.transposed);
            }
            
            assert( maxIterations != -1 );
            
            S = Q.transposed * this;
        }
    
        /**
        Standard operators that have intuitive meaning, same as in classical math.
        
        Note that division operators do no cheks of value of k, so in case of division
        by 0 result matrix will have infinity components. You can check this with isNormal()
        method.
        */
        const pure nothrow @safe Matrix33 opNeg()
        {
            return Matrix33(-m00, -m01, -m02,
                            -m10, -m11, -m12,
                            -m20, -m21, -m22);
        }
    
        /** ditto */
        const pure nothrow @safe Matrix33 opAdd(Matrix33 mat)
        {
            return Matrix33(m00 + mat.m00, m01 + mat.m01, m02 + mat.m02,
                            m10 + mat.m10, m11 + mat.m11, m12 + mat.m12,
                            m20 + mat.m20, m21 + mat.m21, m22 + mat.m22);
        }
    
        /** ditto */
        pure nothrow @safe void opAddAssign(Matrix33 mat)
        {
            m00 += mat.m00; m01 += mat.m01; m02 += mat.m02;
            m10 += mat.m10; m11 += mat.m11; m12 += mat.m12;
            m20 += mat.m20; m21 += mat.m21; m22 += mat.m22;
        }
    
        /** ditto */
        const pure nothrow @safe Matrix33 opSub(Matrix33 mat)
        {
            return Matrix33(m00 - mat.m00, m01 - mat.m01, m02 - mat.m02,
                            m10 - mat.m10, m11 - mat.m11, m12 - mat.m12,
                            m20 - mat.m20, m21 - mat.m21, m22 - mat.m22);
        }
    
        /** ditto */
        pure nothrow @safe void opSubAssign(Matrix33 mat)
        {
            m00 -= mat.m00; m01 -= mat.m01; m02 -= mat.m02;
            m10 -= mat.m10; m11 -= mat.m11; m12 -= mat.m12;
            m20 -= mat.m20; m21 -= mat.m21; m22 -= mat.m22;        
        }
    
        /** ditto */
        const pure nothrow @safe Matrix33 opMul(float_t k)
        {
            return Matrix33(m00 * k, m01 * k, m02 * k,
                            m10 * k, m11 * k, m12 * k,
                            m20 * k, m21 * k, m22 * k);
        }
    
        /** ditto */
        pure nothrow @safe void opMulAssign(float_t k)
        {
            m00 *= k; m01 *= k; m02 *= k;
            m10 *= k; m11 *= k; m12 *= k;
            m20 *= k; m21 *= k; m22 *= k;
        }
    
        /** ditto */
        const pure nothrow @safe Matrix33 opMul_r(float_t k)
        {
            return Matrix33(m00 * k, m01 * k, m02 * k,
                            m10 * k, m11 * k, m12 * k,
                            m20 * k, m21 * k, m22 * k);
        }
    
        /** ditto */
        const pure nothrow @safe Matrix33 opDiv(float_t k)
        {
            
            return Matrix33(m00 / k, m01 / k, m02 / k,
                            m10 / k, m11 / k, m12 / k,
                            m20 / k, m21 / k, m22 / k);
        }
    
        /** ditto */
        pure nothrow @safe void opDivAssign(float_t k)
        {
            m00 /= k; m01 /= k; m02 /= k;
            m10 /= k; m11 /= k; m12 /= k;
            m20 /= k; m21 /= k; m22 /= k;
        }
    
        /** ditto */
        const pure nothrow @safe bool opEquals(Matrix33 mat)
        {
            return m00 == mat.m00 && m01 == mat.m01 && m02 == mat.m02 &&
                   m10 == mat.m10 && m11 == mat.m11 && m12 == mat.m12 &&
                   m20 == mat.m20 && m21 == mat.m21 && m22 == mat.m22;
        }

        /** ditto */
        const pure nothrow @safe Matrix33 opMul(Matrix33 mat)
        {
            return Matrix33(m00 * mat.m00 + m01 * mat.m10 + m02 * mat.m20,
                            m00 * mat.m01 + m01 * mat.m11 + m02 * mat.m21,
                            m00 * mat.m02 + m01 * mat.m12 + m02 * mat.m22,
                            m10 * mat.m00 + m11 * mat.m10 + m12 * mat.m20,
                            m10 * mat.m01 + m11 * mat.m11 + m12 * mat.m21,
                            m10 * mat.m02 + m11 * mat.m12 + m12 * mat.m22,
                            m20 * mat.m00 + m21 * mat.m10 + m22 * mat.m20,
                            m20 * mat.m01 + m21 * mat.m11 + m22 * mat.m21,
                            m20 * mat.m02 + m21 * mat.m12 + m22 * mat.m22 );
        }
    
        /** ditto */
        pure nothrow @safe void opMulAssign(Matrix33 mat)
        {
            this = this * mat;
        }
    
        /** ditto */
        const pure nothrow @safe Vector3 opMul(Vector3 v)
        {
            return Vector3(v.x * m00 + v.y * m01 + v.z * m02,
                           v.x * m10 + v.y * m11 + v.z * m12,
                           v.x * m20 + v.y * m21 + v.z * m22 );
        }
    
        /** Returns: Element at row'th _row and col'th column. */
        const pure @safe float_t opIndex(uint row, uint col)
        in { assert( row < 3 && col < 3 ); }
        body
        {
            return m[col][row];
        }
    
        /** Assigns value f to element at row'th _row and col'th column. */
        pure @safe void opIndexAssign(float_t f, uint row, uint col)
        in { assert( row < 3 && col < 3 ); }
        body
        {
            m[col][row] = f;
        }
        
        /** Returns: Vector representing col'th column. */
        const pure @safe Vector3 opIndex(uint col)
        in { assert( col < 3 ); }
        body
        {
            return v[col];
        }
        
        /** Replaces elements in col'th column with v's values. */
        pure @safe Vector3 opIndexAssign(Vector3 v, uint col)
        in { assert( col < 3 ); }
        body
        {
            return this.v[col] = v;
        }
    
        /**
        Returns: float_t pointer to [0,0] element of this matrix. It's like a _ptr method for arrays.
        
        Remember about column-major matrix memory layout.
        */
        pure nothrow @safe float_t* ptr()
        {
            return a.ptr;
        }
        
        /** Returns: Copy of this matrix with float type elements. */
        const pure nothrow @safe Matrix33f toMatrix33f()
        {
            return Matrix33f(
                cast(float)m00, cast(float)m01, cast(float)m02,
                cast(float)m10, cast(float)m11, cast(float)m12,
                cast(float)m20, cast(float)m21, cast(float)m22 );
        }
        
        /** Returns: Copy of this matrix with double type elements. */
        const pure nothrow @safe Matrix33d toMatrix33d()
        {
            return Matrix33d(
                cast(double)m00, cast(double)m01, cast(double)m02,
                cast(double)m10, cast(double)m11, cast(double)m12,
                cast(double)m20, cast(double)m21, cast(double)m22 );
        }
        
        /** Returns: Copy of this matrix with real type elements. */
        const pure nothrow @safe Matrix33r toMatrix33r()
        {
            return Matrix33r(
                cast(real)m00, cast(real)m01, cast(real)m02,
                cast(real)m10, cast(real)m11, cast(real)m12,
                cast(real)m20, cast(real)m21, cast(real)m22 );
        }

        const string toString() { 
            return format("[" ,m00, ", " ,m01, ", " ,m02, ",\n",
                          " " ,m10, ", " ,m11, ", " ,m12, ",\n",
                          " " ,m20, ", " ,m21, ", " ,m22, "]");
        }
    }
    
    
    alias EqualityByNorm!(Matrix33).equal equal; /// Introduces approximate equality function for Matrix33.
    alias Lerp!(Matrix33).lerp lerp;             /// Introduces linear interpolation function for Matrix33.
    
    /************************************************************************************
    4x4 Matrix.

    Helix matrices uses column-major memory layout.
    *************************************************************************************/
    struct Matrix44
    {
        align (1) union
        {
            struct
            {
                float_t m00, m10, m20, m30;
                float_t m01, m11, m21, m31;
                float_t m02, m12, m22, m32;
                float_t m03, m13, m23, m33;
            }
    
            float_t[4][4] m;
            float_t[16]   a;
            Vector4[4]    v;
        }
    
        /// Identity matrix.
        static immutable Matrix44 identity = {
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1 };
        /// Matrix with all elements set to NaN.
        static immutable Matrix44 nan = {
            float_t.nan, float_t.nan, float_t.nan, float_t.nan,
            float_t.nan, float_t.nan, float_t.nan, float_t.nan,
            float_t.nan, float_t.nan, float_t.nan, float_t.nan,
            float_t.nan, float_t.nan, float_t.nan, float_t.nan };
        /// Matrix with all elements set to 0.
        static immutable Matrix44 zero = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0 };
    
        /**
        Methods to construct matrix in C-like syntax.

        In case with array remember about column-major matrix memory layout,
        note last line with assert in example - it'll be passed.

        Examples:
        ------------
        Matrix33 mat1 = Matrix33(
             1,  2,  3,  4,
             5,  6,  7,  8,
             9, 10, 11, 12,
            13, 14, 15, 16 );
            
        static float[16] a = [
             1,  2,  3,  4,
             5,  6,  7,  8,
             9, 10, 11, 12,
            13, 14, 15, 16 ];
        Matrix33 mat2 = Matrix33(a);

        assert(mat1 == mat2.transposed);
        ------------
        */
        pure nothrow @safe static Matrix44 opCall(float_t m00, float_t m01, float_t m02, float_t m03,
                               float_t m10, float_t m11, float_t m12, float_t m13,
                               float_t m20, float_t m21, float_t m22, float_t m23,
                               float_t m30, float_t m31, float_t m32, float_t m33)
        {
            Matrix44 mat;
            mat.m00 = m00;        mat.m01 = m01;        mat.m02 = m02;        mat.m03 = m03;
            mat.m10 = m10;        mat.m11 = m11;        mat.m12 = m12;        mat.m13 = m13;
            mat.m20 = m20;        mat.m21 = m21;        mat.m22 = m22;        mat.m23 = m23;
            mat.m30 = m30;        mat.m31 = m31;        mat.m32 = m32;        mat.m33 = m33;
            return mat;
        }
        
        /** ditto */
        pure @safe static Matrix44 opCall(float_t[16] a)
        {
            Matrix44 mat;
            mat.a[0..16] = a[0..16].dup;
            return mat;
        }
        
        /**
        Method to construct matrix in C-like syntax. Sets columns to passed vector
        arguments.
        */
        pure nothrow @safe static Matrix44 opCall(Vector4 basisX, Vector4 basisY, Vector4 basisZ,
                               Vector4 basisW = Vector4(0, 0, 0, 1))
        {
            Matrix44 mat;
            mat.v[0] = basisX;
            mat.v[1] = basisY;
            mat.v[2] = basisZ;
            mat.v[3] = basisW;
            return mat;
        }
        
        /**
        Method to construct matrix in C-like syntax. Constructs affine transform
        matrix based on passed vector arguments.
        References:
            $(LINK http://en.wikipedia.org/wiki/Affine_transformation).
        */
        pure nothrow @safe static Matrix44 opCall(Vector3 basisX, Vector3 basisY, Vector3 basisZ,
                               Vector3 translation = Vector3(0, 0, 0))
        {
            return opCall(Vector4(basisX, 0), Vector4(basisX, 0), Vector4(basisX, 0), Vector4(translation, 1));
        }
    
        /** Sets elements to passed values. */
        pure nothrow @safe void set(float_t m00, float_t m01, float_t m02, float_t m03,
                 float_t m10, float_t m11, float_t m12, float_t m13,
                 float_t m20, float_t m21, float_t m22, float_t m23,
                 float_t m30, float_t m31, float_t m32, float_t m33)
        {
            this.m00 = m00;        this.m01 = m01;        this.m02 = m02;        this.m03 = m03;
            this.m10 = m10;        this.m11 = m11;        this.m12 = m12;        this.m13 = m13;
            this.m20 = m20;        this.m21 = m21;        this.m22 = m22;        this.m23 = m23;
            this.m30 = m30;        this.m31 = m31;        this.m32 = m32;        this.m33 = m33;    
        }
        
        /** Sets elements as _a copy of a contents. Remember about column-major matrix memory layout. */
        pure @safe void set(float_t[16] a)
        {
            this.a[0..16] = a[0..16].dup;
        }
    
        /** Sets columns to passed basis vectors. */
        pure nothrow @safe void set(Vector4 basisX, Vector4 basisY, Vector4 basisZ,
                 Vector4 basisW = Vector4(0, 0, 0, 1))
        {
            v[0] = basisX;
            v[1] = basisY;
            v[2] = basisZ;
            v[3] = basisW;
        }

        /** Returns: Whether all components are normalized numbers. */
        const pure nothrow @safe bool isNormal()
        {
            return
                std.math.isNormal(m00) && std.math.isNormal(m01) && std.math.isNormal(m02) && std.math.isNormal(m03) &&
                std.math.isNormal(m10) && std.math.isNormal(m11) && std.math.isNormal(m12) && std.math.isNormal(m13) &&
                std.math.isNormal(m20) && std.math.isNormal(m21) && std.math.isNormal(m22) && std.math.isNormal(m23) &&
                std.math.isNormal(m30) && std.math.isNormal(m31) && std.math.isNormal(m32) && std.math.isNormal(m33);
        }

        /**
        Returns: Whether this matrix is identity.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                               Have the same meaning as in equal function.
        */
        const pure nothrow @safe bool isIdentity(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal(this, identity, relprec, absprec);
        }
        
        /**
        Returns: Whether this matrix is zero.
        Params:
            relprec, absprec = Parameters passed to equal function while calculations.
                        Has the same meaning as in equal function.
        */
        const pure nothrow @safe bool isZero(int relprec = defrelprec, int absprec = defabsprec)
        {
            return equal(normSquare(), 0, relprec, absprec);
        }
        
        /**
        Resets this matrix to affine transform matrix based on passed
        vector arguments.
        */
        pure nothrow @safe void set(Vector3 basisX, Vector3 basisY, Vector3 basisZ,
                 Vector3 translation = Vector3(0, 0, 0))
        {
            v[0] = Vector4(basisX, 0);
            v[1] = Vector4(basisY, 0);
            v[2] = Vector4(basisZ, 0);
            v[3] = Vector4(translation, 1);
        }
    
        /** Constructs _scale matrix with _scale coefficients specified as arguments. */
        pure nothrow @safe static Matrix44 scale(float_t x, float_t y, float_t z)
        {
            Matrix44 mat = identity;
            with (mat)
            {
                m00 = x;
                m11 = y;
                m22 = z;            
            }
    
            return mat;
        }
    
        /** Constructs _scale matrix with _scale coefficients specified as v's components. */
        pure nothrow @safe static Matrix44 scale(Vector3 v)
        {
            return scale(v.x, v.y, v.z);
        }
    
        /** Construct matrix that represents rotation around corresponding axis. */
        pure nothrow @safe static Matrix44 rotationX(float_t radians)
        {
            Matrix44 mat = identity;
            float_t c = cos(radians);
            float_t s = sin(radians);
            with (mat)
            {
                m11 = m22 = c;
                m21 = s;
                m12 = -s;            
            }
    
            return mat;
        }
    
        /** ditto */
        pure nothrow @safe static Matrix44 rotationY(float_t radians)
        {
            Matrix44 mat = identity;
            float_t c = cos(radians);
            float_t s = sin(radians);
            with (mat)
            {
                m00 = m22 = c;
                m20 = -s;
                m02 = s;            
            }
    
            return mat;
        }
    
        /** ditto */
        pure nothrow @safe static Matrix44 rotationZ(float_t radians)
        {
            Matrix44 mat = identity;
            float_t c = cos(radians);
            float_t s = sin(radians);
            with (mat)
            {
                m00 = m11 = c;
                m10 = s;
                m01 = -s;            
            }
    
            return mat;
        }
    
        /**
        Constructs matrix that represents _rotation specified by euler angles passed as arguments.
        Order of _rotation application is: roll (Z axis), pitch (X axis), yaw (Y axis).
        */
        pure nothrow @safe static Matrix44 rotation(float_t yaw, float_t pitch, float_t roll)
        {
            return Matrix44.rotationY(yaw) * Matrix44.rotationX(pitch) * Matrix44.rotationZ(roll);
        }
    
        /**
        Constructs matrix that represents _rotation specified by axis and angle.
        Method works with assumption that axis is unit vector.        
        Throws:
            AssertError on non-unit axis call attempt if module was compiled with
            contract checks enabled.
        */
        pure @safe static Matrix44 rotation(Vector3 axis, float_t radians)
        in { assert( axis.isUnit() ); }
        body
        {
            real c = cos(radians);
            real s = sin(radians);
            real cc = 1.0 - c;
            real x2 = axis.x * axis.x;
            real y2 = axis.y * axis.y;
            real z2 = axis.z * axis.z;
            real xycc = axis.x * axis.y * cc;
            real xzcc = axis.x * axis.z * cc;
            real yzcc = axis.y * axis.z * cc;
            real xs = axis.x * s;
            real ys = axis.y * s;
            real zs = axis.z * s;
    
            Matrix44 mat = identity;
            with (mat)
            {
                m00 = x2 * cc + c;      m01 = xycc - zs;        m02 = xzcc + ys;
                m10 = xycc + zs;        m11 = y2 * cc + c;      m12 = yzcc - xs;
                m20 = xzcc - ys;        m21 = yzcc + xs;        m22 = z2 * cc + c;
            }
    
            return mat;
        }
        
        /**
        Constructs matrix that represents _rotation specified by quaternion.
        Method works with assumption that quaternion is unit.        
        Throws:
            AssertError on non-unit quaternion call attempt if module was compiled with
            contract checks enabled.
        */
        pure @safe static Matrix44 rotation(Quaternion q)
        in { assert( q.isUnit() ); }
        body
        {
            float_t tx  = 2.0f * q.x;
            float_t ty  = 2.0f * q.y;
            float_t tz  = 2.0f * q.z;
            float_t twx = tx * q.w;
            float_t twy = ty * q.w;
            float_t twz = tz * q.w;
            float_t txx = tx * q.x;
            float_t txy = ty * q.x;
            float_t txz = tz * q.x;
            float_t tyy = ty * q.y;
            float_t tyz = tz * q.y;
            float_t tzz = tz * q.z;
            
            Matrix44 mat = identity;
            with (mat)
            {
                m00 = 1.0f - (tyy + tzz); m01 = txy - twz;          m02 = txz + twy;
                m10 = txy + twz;          m11 = 1.0f - (txx + tzz); m12 = tyz - twx;
                m20 = txz - twy;          m21 = tyz + twx;          m22 = 1.0f - (txx + tyy);
            }
            
            return mat;
        }
    
        /** Constructs _translation matrix with offset values specified as arguments. */
        pure nothrow @safe static Matrix44 translation(float_t x, float_t y, float_t z)
        {
            return Matrix44(1, 0, 0, x,
                            0, 1, 0, y,
                            0, 0, 1, z,
                            0, 0, 0, 1);
        }
    
        /** Constructs _translation matrix with offset values specified as v's components. */
        pure nothrow @safe static Matrix44 translation(Vector3 v)
        {
            return translation(v.x, v.y, v.z);
        }
        
        /**
        Constructs one-point perspecive projection matrix.
        Params:
            fov =       Field of view in vertical plane in radians.
            aspect =    Frustum's width / height coefficient. It shouldn't be 0.
            near =      Distance to near plane.
            near =      Distance to far plane.
        */
        pure @safe static Matrix44 perspective(float_t fov, float_t aspect, float_t near, float_t far)
        in
        {
            assert( fov < 2*PI );
            assert( !equal(aspect, 0) );
            assert( near > 0 );
            assert( far > near );
        }
        body
        {
            real cot = 1. / tan(fov / 2.);
                    
            return Matrix44(cot / aspect,    0,                            0,                                  0,
                            0,             cot,                            0,                                  0,
                            0,               0,  (near + far) / (near - far), 2.0f * (near * far) / (near - far),
                            0,               0,                           -1,                                  0);
        }
        
        /**
        Constructs view matrix.
        Params:
            eye =       Viewer's eye position.
            target =    View target.
            up =        View up vector.
        
        Arguments should not be complanar, elsewise matrix will contain infinity
        elements. You can check this with isNormal() method.
        */
        pure nothrow @safe static Matrix44 lookAt(Vector3 eye, Vector3 target, Vector3 up)
        {
            Vector3 z = (eye - target).normalized();
            alias up y;
            Vector3 x = cross(y, z);
            y = cross(z, x);
            x.normalize();
            y.normalize();
                    
            Matrix44 mat = identity;
            mat.v[0].xyz = Vector3(x.x, y.x, z.x);
            mat.v[1].xyz = Vector3(x.y, y.y, z.y);
            mat.v[2].xyz = Vector3(x.z, y.z, z.z);
                    
            mat.m03 = -dot(eye, x);
            mat.m13 = -dot(eye, y);
            mat.m23 = -dot(eye, z);
                    
            return mat;    
        }
        
        /**
        Returns: Inverse copy of this matrix.
        
        In case if this matrix is singular (i.e. determinant = 0) result matrix will has
        infinity elements. You can check this with isNormal() method.
        */
        const pure nothrow @safe @property Matrix44 inverse()
        {
            real det = determinant();
            //if (equal(det, 0))
            //{
            //    return nan;
            //}
            
            real rdet = 1/det;
            return Matrix44(
                rdet * (m11 * (m22 * m33 - m23 * m32) + m12 * (m23 * m31 - m21 * m33) + m13 * (m21 * m32 - m22 * m31)),
                rdet * (m21 * (m02 * m33 - m03 * m32) + m22 * (m03 * m31 - m01 * m33) + m23 * (m01 * m32 - m02 * m31)),
                rdet * (m31 * (m02 * m13 - m03 * m12) + m32 * (m03 * m11 - m01 * m13) + m33 * (m01 * m12 - m02 * m11)),
                rdet * (m01 * (m13 * m22 - m12 * m23) + m02 * (m11 * m23 - m13 * m21) + m03 * (m12 * m21 - m11 * m22)),
                rdet * (m12 * (m20 * m33 - m23 * m30) + m13 * (m22 * m30 - m20 * m32) + m10 * (m23 * m32 - m22 * m33)),
                rdet * (m22 * (m00 * m33 - m03 * m30) + m23 * (m02 * m30 - m00 * m32) + m20 * (m03 * m32 - m02 * m33)),
                rdet * (m32 * (m00 * m13 - m03 * m10) + m33 * (m02 * m10 - m00 * m12) + m30 * (m03 * m12 - m02 * m13)),
                rdet * (m02 * (m13 * m20 - m10 * m23) + m03 * (m10 * m22 - m12 * m20) + m00 * (m12 * m23 - m13 * m22)),
                rdet * (m13 * (m20 * m31 - m21 * m30) + m10 * (m21 * m33 - m23 * m31) + m11 * (m23 * m30 - m20 * m33)),
                rdet * (m23 * (m00 * m31 - m01 * m30) + m20 * (m01 * m33 - m03 * m31) + m21 * (m03 * m30 - m00 * m33)),
                rdet * (m33 * (m00 * m11 - m01 * m10) + m30 * (m01 * m13 - m03 * m11) + m31 * (m03 * m10 - m00 * m13)),
                rdet * (m03 * (m11 * m20 - m10 * m21) + m00 * (m13 * m21 - m11 * m23) + m01 * (m10 * m23 - m13 * m20)),
                rdet * (m10 * (m22 * m31 - m21 * m32) + m11 * (m20 * m32 - m22 * m30) + m12 * (m21 * m30 - m20 * m31)),
                rdet * (m20 * (m02 * m31 - m01 * m32) + m21 * (m00 * m32 - m02 * m30) + m22 * (m01 * m30 - m00 * m31)),
                rdet * (m30 * (m02 * m11 - m01 * m12) + m31 * (m00 * m12 - m02 * m10) + m32 * (m01 * m10 - m00 * m11)),
                rdet * (m00 * (m11 * m22 - m12 * m21) + m01 * (m12 * m20 - m10 * m22) + m02 * (m10 * m21 - m11 * m20)));
        }
        
        /**
        Inverts this matrix.
        
        In case if matrix is singular (i.e. determinant = 0) result matrix will has
        infinity elements. You can check this with isNormal() method.
        */
        pure nothrow @safe void invert()
        {
            real det = determinant();
            //if (equal(det, 0))
            //{
            //    *this = nan;
            //    return;
            //}
            
            real rdet = 1/det;
            set(rdet * (m11 * (m22 * m33 - m23 * m32) + m12 * (m23 * m31 - m21 * m33) + m13 * (m21 * m32 - m22 * m31)),
                rdet * (m21 * (m02 * m33 - m03 * m32) + m22 * (m03 * m31 - m01 * m33) + m23 * (m01 * m32 - m02 * m31)),
                rdet * (m31 * (m02 * m13 - m03 * m12) + m32 * (m03 * m11 - m01 * m13) + m33 * (m01 * m12 - m02 * m11)),
                rdet * (m01 * (m13 * m22 - m12 * m23) + m02 * (m11 * m23 - m13 * m21) + m03 * (m12 * m21 - m11 * m22)),
                rdet * (m12 * (m20 * m33 - m23 * m30) + m13 * (m22 * m30 - m20 * m32) + m10 * (m23 * m32 - m22 * m33)),
                rdet * (m22 * (m00 * m33 - m03 * m30) + m23 * (m02 * m30 - m00 * m32) + m20 * (m03 * m32 - m02 * m33)),
                rdet * (m32 * (m00 * m13 - m03 * m10) + m33 * (m02 * m10 - m00 * m12) + m30 * (m03 * m12 - m02 * m13)),
                rdet * (m02 * (m13 * m20 - m10 * m23) + m03 * (m10 * m22 - m12 * m20) + m00 * (m12 * m23 - m13 * m22)),
                rdet * (m13 * (m20 * m31 - m21 * m30) + m10 * (m21 * m33 - m23 * m31) + m11 * (m23 * m30 - m20 * m33)),
                rdet * (m23 * (m00 * m31 - m01 * m30) + m20 * (m01 * m33 - m03 * m31) + m21 * (m03 * m30 - m00 * m33)),
                rdet * (m33 * (m00 * m11 - m01 * m10) + m30 * (m01 * m13 - m03 * m11) + m31 * (m03 * m10 - m00 * m13)),
                rdet * (m03 * (m11 * m20 - m10 * m21) + m00 * (m13 * m21 - m11 * m23) + m01 * (m10 * m23 - m13 * m20)),
                rdet * (m10 * (m22 * m31 - m21 * m32) + m11 * (m20 * m32 - m22 * m30) + m12 * (m21 * m30 - m20 * m31)),
                rdet * (m20 * (m02 * m31 - m01 * m32) + m21 * (m00 * m32 - m02 * m30) + m22 * (m01 * m30 - m00 * m31)),
                rdet * (m30 * (m02 * m11 - m01 * m12) + m31 * (m00 * m12 - m02 * m10) + m32 * (m01 * m10 - m00 * m11)),
                rdet * (m00 * (m11 * m22 - m12 * m21) + m01 * (m12 * m20 - m10 * m22) + m02 * (m10 * m21 - m11 * m20)));
        }
        
        /** Returns: Determinant */
        const pure nothrow @safe @property real determinant()
        {
            return
                + (m00 * m11 - m01 * m10) * (m22 * m33 - m23 * m32)
                - (m00 * m12 - m02 * m10) * (m21 * m33 - m23 * m31)
                + (m00 * m13 - m03 * m10) * (m21 * m32 - m22 * m31)
                + (m01 * m12 - m02 * m11) * (m20 * m33 - m23 * m30)
                - (m01 * m13 - m03 * m11) * (m20 * m32 - m22 * m30)
                + (m02 * m13 - m03 * m12) * (m20 * m31 - m21 * m30);
        }
        
        /**
        Returns: Frobenius _norm of matrix.
        References:
            $(LINK http://en.wikipedia.org/wiki/Frobenius_norm#Frobenius_norm).        
        */
        const pure nothrow @safe @property real norm()
        {
            return sqrt( normSquare );
        }
        
        /**
        Returns: Square of Frobenius norm of matrix.

        Method doesn't need calculation of square root.

        References:
            $(LINK http://en.wikipedia.org/wiki/Frobenius_norm#Frobenius_norm).
        */
        const pure nothrow @safe @property real normSquare()
        {
            real ret = 0;
            for (int i = 16; i--; )
            {
                real x = a[i];
                ret += x * x;
            }
            
            return ret;
        }
        
        /** 
        Returns: Whether this matrix represents affine transformation.
        References:
            $(LINK http://en.wikipedia.org/wiki/Affine_transformation).
        */
        const pure nothrow @safe bool isAffine()
        {
            return equal(m30, 0) && equal(m31, 0) && equal(m32, 0) && equal(m33, 1);
        }
        
        /** Transposes this matrix. */
        pure nothrow @safe void transpose()
        {
            /*           */        swap(m01, m10);        swap(m02, m20);        swap(m03, m30);
            /*           */        /*           */        swap(m12, m21);        swap(m13, m31);
            /*           */        /*           */        /*           */        swap(m23, m32);
            /*           */        /*           */        /*           */        /*           */
        }
        
        /** Returns: Transposed copy of this matrix. */
        const pure nothrow @safe @property Matrix44 transposed()
        {
            return Matrix44(
                m00, m10, m20, m30,
                m01, m11, m21, m31,
                m02, m12, m22, m32,
                m03, m13, m23, m33 );
        }
        
        /** R/W property. Corner 3x3 minor. */
        const pure nothrow @safe @property Matrix33 cornerMinor()
        {
            return Matrix33(m00, m01, m02,
                            m10, m11, m12,
                            m20, m21, m22);
        }
        
        /** ditto */
        pure nothrow @safe @property void cornerMinor(Matrix33 mat)
        {
            m00 = mat.m00;        m01 = mat.m01;        m02 = mat.m02;
            m10 = mat.m10;        m11 = mat.m11;        m12 = mat.m12;
            m20 = mat.m20;        m21 = mat.m21;        m22 = mat.m22;
        }
        
        /**
        Standard operators that have intuitive meaning, same as in classical math. Exception
        is multiplication with Vecto3 that doesn't make sense for classical math, in that case
        Vector3 is implicitl expanded to Vector4 with w=1.
        
        Note that division operators do no cheks of value of k, so in case of division
        by 0 result matrix will have infinity components. You can check this with isNormal()
        method.
        */
        const pure nothrow @safe Matrix44 opNeg()
        {
            return Matrix44(-m00, -m01, -m02, -m03,
                            -m10, -m11, -m12, -m13,
                            -m20, -m21, -m22, -m23,
                            -m30, -m31, -m32, -m33);
        }
    
        /** ditto */
        const pure nothrow @safe Matrix44 opAdd(Matrix44 mat)
        {
            return Matrix44(m00 + mat.m00, m01 + mat.m01, m02 + mat.m02, m03 + mat.m03,
                            m10 + mat.m10, m11 + mat.m11, m12 + mat.m12, m13 + mat.m13,
                            m20 + mat.m20, m21 + mat.m21, m22 + mat.m22, m23 + mat.m23,
                            m30 + mat.m30, m31 + mat.m31, m32 + mat.m32, m33 + mat.m33);
        }
    
        /** ditto */
        pure nothrow @safe void opAddAssign(Matrix44 mat)
        {
            m00 += mat.m00; m01 += mat.m01; m02 += mat.m02; m03 += mat.m03;
            m10 += mat.m10; m11 += mat.m11; m12 += mat.m12; m13 += mat.m13;
            m20 += mat.m20; m21 += mat.m21; m22 += mat.m22; m23 += mat.m23;
            m30 += mat.m30; m31 += mat.m31; m32 += mat.m32; m33 += mat.m33;
        }
    
        /** ditto */
        const pure nothrow @safe Matrix44 opSub(Matrix44 mat)
        {
            return Matrix44(m00 - mat.m00, m01 - mat.m01, m02 - mat.m02, m03 - mat.m03,
                            m10 - mat.m10, m11 - mat.m11, m12 - mat.m12, m13 - mat.m13,
                            m20 - mat.m20, m21 - mat.m21, m22 - mat.m22, m23 - mat.m23,
                            m30 - mat.m30, m31 - mat.m31, m32 - mat.m32, m33 - mat.m33);
        }
    
        /** ditto */
        pure nothrow @safe void opSubAssign(Matrix44 mat)
        {
            m00 -= mat.m00; m01 -= mat.m01; m02 -= mat.m02; m03 -= mat.m03;
            m10 -= mat.m10; m11 -= mat.m11; m12 -= mat.m12; m13 -= mat.m13;
            m20 -= mat.m20; m21 -= mat.m21; m22 -= mat.m22; m23 -= mat.m23;        
            m30 -= mat.m30; m31 -= mat.m31; m32 -= mat.m32; m33 -= mat.m33;        
        }
    
        /** ditto */
        const pure nothrow @safe Matrix44 opMul(float_t k)
        {
            return Matrix44(m00 * k, m01 * k, m02 * k, m03 * k,
                            m10 * k, m11 * k, m12 * k, m13 * k,
                            m20 * k, m21 * k, m22 * k, m23 * k,
                            m30 * k, m31 * k, m32 * k, m33 * k);
        }
    
        /** ditto */
        pure nothrow @safe void opMulAssign(float_t k)
        {
            m00 *= k; m01 *= k; m02 *= k; m03 *= k;
            m10 *= k; m11 *= k; m12 *= k; m13 *= k;
            m20 *= k; m21 *= k; m22 *= k; m23 *= k;
            m30 *= k; m31 *= k; m32 *= k; m33 *= k;
        }
    
        /** ditto */
        const pure nothrow @safe Matrix44 opMul_r(float_t k)
        {
            return Matrix44(m00 * k, m01 * k, m02 * k, m03 * k,
                            m10 * k, m11 * k, m12 * k, m13 * k,
                            m20 * k, m21 * k, m22 * k, m23 * k,
                            m30 * k, m31 * k, m32 * k, m33 * k);
        }
    
        /** ditto */
        const pure nothrow @safe Matrix44 opDiv(float_t k)
        {
            
            return Matrix44(m00 / k, m01 / k, m02 / k, m03 / k,
                            m10 / k, m11 / k, m12 / k, m13 / k,
                            m20 / k, m21 / k, m22 / k, m23 / k,
                            m30 / k, m31 / k, m32 / k, m33 / k);
        }
    
        /** ditto */
        pure nothrow @safe void opDivAssign(float_t k)
        {
            m00 /= k; m01 /= k; m02 /= k; m03 /= k;
            m10 /= k; m11 /= k; m12 /= k; m13 /= k;
            m20 /= k; m21 /= k; m22 /= k; m23 /= k;
            m30 /= k; m31 /= k; m32 /= k; m33 /= k;
        }
    
        /** ditto */
        const pure nothrow @safe bool opEquals(Matrix44 mat)
        {
            return m00 == mat.m00 && m01 == mat.m01 && m02 == mat.m02 && m03 == mat.m03 &&
                   m10 == mat.m10 && m11 == mat.m11 && m12 == mat.m12 && m13 == mat.m13 &&
                   m20 == mat.m20 && m21 == mat.m21 && m22 == mat.m22 && m23 == mat.m23 &&
                   m30 == mat.m30 && m31 == mat.m31 && m32 == mat.m32 && m33 == mat.m33;
        }

        /** ditto */
        const pure nothrow @safe Matrix44 opMul(Matrix44 mat)
        {
            return Matrix44(m00 * mat.m00 + m01 * mat.m10 + m02 * mat.m20 + m03 * mat.m30,
                            m00 * mat.m01 + m01 * mat.m11 + m02 * mat.m21 + m03 * mat.m31,
                            m00 * mat.m02 + m01 * mat.m12 + m02 * mat.m22 + m03 * mat.m32,
                            m00 * mat.m03 + m01 * mat.m13 + m02 * mat.m23 + m03 * mat.m33,
    
                            m10 * mat.m00 + m11 * mat.m10 + m12 * mat.m20 + m13 * mat.m30,
                            m10 * mat.m01 + m11 * mat.m11 + m12 * mat.m21 + m13 * mat.m31,
                            m10 * mat.m02 + m11 * mat.m12 + m12 * mat.m22 + m13 * mat.m32,
                            m10 * mat.m03 + m11 * mat.m13 + m12 * mat.m23 + m13 * mat.m33,
    
                            m20 * mat.m00 + m21 * mat.m10 + m22 * mat.m20 + m23 * mat.m30,
                            m20 * mat.m01 + m21 * mat.m11 + m22 * mat.m21 + m23 * mat.m31,
                            m20 * mat.m02 + m21 * mat.m12 + m22 * mat.m22 + m23 * mat.m32,
                            m20 * mat.m03 + m21 * mat.m13 + m22 * mat.m23 + m23 * mat.m33,
    
                            m30 * mat.m00 + m31 * mat.m10 + m32 * mat.m20 + m33 * mat.m30,
                            m30 * mat.m01 + m31 * mat.m11 + m32 * mat.m21 + m33 * mat.m31,
                            m30 * mat.m02 + m31 * mat.m12 + m32 * mat.m22 + m33 * mat.m32,
                            m30 * mat.m03 + m31 * mat.m13 + m32 * mat.m23 + m33 * mat.m33);
        }
    
        /** ditto */
        pure nothrow @safe void opMulAssign(Matrix44 mat)
        {
            this = this * mat;
        }
    
        /** ditto */
        const pure nothrow @safe Vector3 opMul(Vector3 v)
        {
            return Vector3(v.x * m00 + v.y * m01 + v.z * m02 + m03,
                           v.x * m10 + v.y * m11 + v.z * m12 + m13,
                           v.x * m20 + v.y * m21 + v.z * m22 + m23 );
        }
    
        /** ditto */
        const pure nothrow @safe Vector4 opMul(Vector4 v)
        {
            return Vector4(v.x * m00 + v.y * m01 + v.z * m02 + v.w * m03,
                           v.x * m10 + v.y * m11 + v.z * m12 + v.w * m13,
                           v.x * m20 + v.y * m21 + v.z * m22 + v.w * m23,
                           v.x * m30 + v.y * m31 + v.z * m32 + v.w * m33);
        }
    
        /** Returns: Element at row'th _row and col'th column. */
        const pure @safe float_t opIndex(uint row, uint col)
        in { assert( col < 4 && row < 4 ); }
        body
        {
            return m[col][row];
        }
    
        /** Assigns value f to element at row'th _row and col'th column. */
        pure @safe void opIndexAssign(float_t f, uint row, uint col)
        in { assert( col < 4 && row < 4 ); }
        body
        {
            m[col][row] = f;
        }
        
        /** Returns: Vector representing col'th column. */
        const pure @safe Vector4 opIndex(uint col)
        in { assert( col < 4 ); }
        body
        {
            return v[col];
        }
    
        /** Replaces elements in col'th column with v's values. */
        pure @safe void opIndexAssign(Vector4 v, uint col)
        in { assert( col < 4 ); }
        body
        {
            this.v[col] = v;
        }
    
        /**
        Returns: float_t pointer to [0,0] element of this matrix. It's like a _ptr method for arrays.
        
        Remember about column-major matrix memory layout.
        */
        pure @safe float_t* ptr()
        {
            return a.ptr;
        }
        
        /** Returns: Copy of this matrix with float type elements. */
        const pure nothrow @safe Matrix44f toMatrix44f()
        {
            return Matrix44f(
                cast(float)m00, cast(float)m01, cast(float)m02, cast(float)m03,
                cast(float)m10, cast(float)m11, cast(float)m12, cast(float)m13,
                cast(float)m20, cast(float)m21, cast(float)m22, cast(float)m23,
                cast(float)m30, cast(float)m31, cast(float)m32, cast(float)m33 );
        }
        
        /** Returns: Copy of this matrix with double type elements. */
        const pure nothrow @safe Matrix44d toMatrix44d()
        {
            return Matrix44d(
                cast(double)m00, cast(double)m01, cast(double)m02, cast(double)m03,
                cast(double)m10, cast(double)m11, cast(double)m12, cast(double)m13,
                cast(double)m20, cast(double)m21, cast(double)m22, cast(double)m23,
                cast(double)m30, cast(double)m31, cast(double)m32, cast(double)m33 );
        }
        
        /** Returns: Copy of this matrix with real type elements. */
        const pure nothrow @safe Matrix44r toMatrix44r()
        {
            return Matrix44r(
                cast(real)m00, cast(real)m01, cast(real)m02, cast(real)m03,
                cast(real)m10, cast(real)m11, cast(real)m12, cast(real)m13,
                cast(real)m20, cast(real)m21, cast(real)m22, cast(real)m23,
                cast(real)m30, cast(real)m31, cast(real)m32, cast(real)m33 );
        }

        const string toString() { 
            return format("[" ,m00, ", " ,m01, ", " ,m02, ", " ,m03, ",\n",
                          " " ,m10, ", " ,m11, ", " ,m12, ", " ,m13, ",\n",
                          " " ,m20, ", " ,m21, ", " ,m22, ", " ,m23, ",\n",
                          " " ,m30, ", " ,m31, ", " ,m32, ", " ,m33, "]");
        }

    }
    
    alias EqualityByNorm!(Matrix44).equal equal; /// Introduces approximate equality function for Matrix44.
    alias Lerp!(Matrix44).lerp lerp;             /// Introduces linear interpolation function for Matrix44.    
}

alias LinearAlgebra!(float).Vector2             Vector2f;
alias LinearAlgebra!(float).Vector3             Vector3f;
alias LinearAlgebra!(float).Vector4             Vector4f;
alias LinearAlgebra!(float).Quaternion          Quaternionf;
alias LinearAlgebra!(float).Matrix22            Matrix22f;
alias LinearAlgebra!(float).Matrix33            Matrix33f;
alias LinearAlgebra!(float).Matrix44            Matrix44f;
alias LinearAlgebra!(float).equal               equal;
alias LinearAlgebra!(float).dot                 dot;
public alias LinearAlgebra!(float).outer        outer;
alias LinearAlgebra!(float).cross               cross;
alias LinearAlgebra!(float).isBasisOrthogonal   isBasisOrthogonal;
alias LinearAlgebra!(float).isBasisOrthonormal  isBasisOrthonormal;
alias LinearAlgebra!(float).lerp                lerp;
alias LinearAlgebra!(float).slerp               slerp;

alias LinearAlgebra!(double).Vector2            Vector2d;
alias LinearAlgebra!(double).Vector3            Vector3d;
alias LinearAlgebra!(double).Vector4            Vector4d;
alias LinearAlgebra!(double).Quaternion         Quaterniond;
alias LinearAlgebra!(double).Matrix22           Matrix22d;
alias LinearAlgebra!(double).Matrix33           Matrix33d;
alias LinearAlgebra!(double).Matrix44           Matrix44d;
alias LinearAlgebra!(double).equal              equal;
alias LinearAlgebra!(double).dot                dot;
//alias LinearAlgebra!(double).outer              outer;
alias LinearAlgebra!(double).cross              cross;
alias LinearAlgebra!(double).isBasisOrthogonal  isBasisOrthogonal;
alias LinearAlgebra!(double).isBasisOrthonormal isBasisOrthonormal;
alias LinearAlgebra!(double).lerp               lerp;
alias LinearAlgebra!(double).slerp              slerp;

alias LinearAlgebra!(real).Vector2              Vector2r;
alias LinearAlgebra!(real).Vector3              Vector3r;
alias LinearAlgebra!(real).Vector4              Vector4r;
alias LinearAlgebra!(real).Quaternion           Quaternionr;
alias LinearAlgebra!(real).Matrix22             Matrix22r;
alias LinearAlgebra!(real).Matrix33             Matrix33r;
alias LinearAlgebra!(real).Matrix44             Matrix44r;
alias LinearAlgebra!(real).equal                equal;
alias LinearAlgebra!(real).dot                  dot;
//alias LinearAlgebra!(real).outer                outer;
alias LinearAlgebra!(real).cross                cross;
alias LinearAlgebra!(real).isBasisOrthogonal    isBasisOrthogonal;
alias LinearAlgebra!(real).isBasisOrthonormal   isBasisOrthonormal;
alias LinearAlgebra!(real).lerp                 lerp;
alias LinearAlgebra!(real).slerp                slerp;

alias LinearAlgebra!(helix.config.float_t).Vector2     Vector2;
alias LinearAlgebra!(helix.config.float_t).Vector3     Vector3;
alias LinearAlgebra!(helix.config.float_t).Vector4     Vector4;
alias LinearAlgebra!(helix.config.float_t).Quaternion  Quaternion;
alias LinearAlgebra!(helix.config.float_t).Matrix22    Matrix22;
alias LinearAlgebra!(helix.config.float_t).Matrix33    Matrix33;
alias LinearAlgebra!(helix.config.float_t).Matrix44    Matrix44;

unittest
{
    assert( Vector2(1, 2).normalized().isUnit() );
    assert( Vector3(1, 2, 3).normalized().isUnit() );
    assert( Vector4(1, 2, 3, 4).normalized().isUnit() );

    assert( Vector2(1, 2).dominatingAxis() == Ort.Y );
    assert( Vector3(1, 2, 3).dominatingAxis() == Ort.Z );
    assert( Vector4(1, 2, 3, 4).dominatingAxis() == Ort.W );

    Vector4 v;
    v.set(1, 2, 3, 4);
    assert( v.isNormal() );
    v /= 0;
    assert( !v.isNormal() );

    v.set(1, 2, 3, 4);
    v[Ort.Y] = v[Ort.X];
    assert( v == Vector4(1, 1, 3, 4) );

    Vector4 t = Vector4(100, 200, 300, 400);
    Vector4 s;
    v.set(1, 2, 3, 4);
    s = v;
    v += t;
    v -= t;
    v = (v + t) - t;
    v *= 100;
    v /= 100;
    v = (10 * v * 10) / 100;
    assert( equal(v, s) );

    assert( dot( cross( Vector3(1, 0, 2), Vector3(4, 0, 5) ), Vector3(3, 0, -2) )  == 0 );
}

unittest
{
    real yaw = PI / 8;
    real pitch = PI / 3;
    real roll = PI / 4;
    
    Quaternion q = Quaternion( Matrix33.rotation(yaw, pitch, roll) );
    assert( equal(q.yaw, yaw) );
    assert( equal(q.pitch, pitch) );
    assert( equal(q.roll, roll) );
}

unittest
{
    Matrix33 mat1 = Matrix33(1,2,3,4,5,6,7,8,9);
    static float[9] a = [1,2,3,4,5,6,7,8,9];
    Matrix33 mat2 = Matrix33(a);

    assert(mat1 == mat2.transposed);
}

/*
unittest
{
    Matrix33 a;
    
    a.m01 = 2;
    a.a[1] = 3;
    a.v[0].z = 4;
    assert(a[0, 1] == 2);
    assert(a[1, 0] == 3);
    assert(a[2, 0] == 4);
}
*/
unittest
{
    Matrix33 a = Matrix33.rotation( Vector3(1, 2, 3).normalized, PI / 7f );
    Matrix33 b = a.inverse;
    b.invert();
    assert( equal(a, b) );
    assert( equal(a.transposed.inverse, a.inverse.transposed) );
}

unittest
{
    Matrix33 Q, S;
    Matrix33 rot = Matrix33.rotationZ(PI / 7);
    Matrix33 scale = Matrix33.scale(-1, 2, 3);
    Matrix33 composition = rot * scale;
    composition.polarDecomposition(Q, S);    
    assert( equal(Q * S, composition) );
}
